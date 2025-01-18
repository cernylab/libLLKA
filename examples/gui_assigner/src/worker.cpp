/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "worker.h"

#include <QEventLoop>
#include <QMessageBox>
#include <QThread>

Worker::Worker(const LLKA::ClassificationContext &m_ctx, std::filesystem::path pathToCif) :
    m_pathToCif{std::move(pathToCif)},
    m_cancel{false},
    m_failed{false},
    m_ctx{m_ctx}
{
}

void Worker::cancel()
{
    m_cancel = true;
}

void Worker::work()
{
    QEventLoop loop{};

    try {
        emit reportProgress(AssignmentInProgressDlg::ReadingCif, 0);

        auto stru = loadStructure(m_pathToCif);
        const auto &stepsToAssign = stru.steps;

        if (stepsToAssign.empty()) {
            m_failed = true;
            emit done({
                WorkerResult::Failure,
                "There are no steps to assign. Did you load the correct file?"
            });
            return;
        }

        std::vector<Step> assignedSteps;

        int percentDone = 0;
        for (size_t idx = 0; idx < stepsToAssign.size(); idx++) {
            const auto &step = stepsToAssign[idx];
            if (m_cancel) {
                emit done({ WorkerResult::Canceled, "" });
                return;
            }

            auto res = LLKA::classifyStep(step.structure, m_ctx);
            if (!res.isSuccess()) {
                assignedSteps.push_back(step);
                assignedSteps.back().state = Step::ClassificationFailed;
            } else {
                const auto &classification = res.success();
                assignedSteps.emplace_back(stru.hasMultipleModels, stru.id, step.structure, classification);
            }

            int perc = 100 * idx / stepsToAssign.size();
            if (perc > percentDone) {
                percentDone = perc;
                emit reportProgress(AssignmentInProgressDlg::Assigning, percentDone);
            }
            loop.processEvents();
        }

        assignedStru = LoadedStructure(stru.hasMultipleModels, std::move(stru.id), std::move(assignedSteps));

        emit done({ WorkerResult::Success, "" });
    } catch (const AssignmentError &err) {
        m_failed = true;
        emit done({
            WorkerResult::Failure,
            QString{"Failed to load structure from Cif file: %1\n%2"}.arg(LLKA::errorToString(err.tRet).c_str()).arg(err.extraMessage.c_str())
        });
    }
}

auto runAssignment(const LLKA::ClassificationContext &ctx, std::filesystem::path pathToCif, QWidget *dialogParent) -> LoadedStructure
{
    auto thread = new QThread{};
    Worker worker{ctx, std::move(pathToCif)};
    worker.moveToThread(thread);

    QString failMessage;
    AssignmentInProgressDlg dlg{dialogParent};

    // Clean up safely when we are done
    QObject::connect(thread, &QThread::finished, thread, &QThread::deleteLater);
    // Begin work when the worker thread starts
    QObject::connect(thread, &QThread::started, &worker, &Worker::work);
    // Wire up dialog to worker signals
    QObject::connect(&worker, &Worker::done, &dlg, &AssignmentInProgressDlg::onDone);
    QObject::connect(&worker, &Worker::reportProgress, &dlg, &AssignmentInProgressDlg::onProgress);
    QObject::connect(&dlg, &AssignmentInProgressDlg::cancel, &worker, &Worker::cancel);

    QObject::connect(&worker, &Worker::done, [thread, &failMessage](const WorkerResult &result) {
        thread->quit();
        if (result.outcome == WorkerResult::Failure)
            failMessage = result.message;
    });

    thread->start();
    dlg.exec();
    thread->wait();

    if (!worker.suceeded()) {
        if (!failMessage.isEmpty()) {
            QMessageBox::warning(
                nullptr,
                "Assigment failed",
                failMessage,
                QMessageBox::Ok
            );
        }
        return {};
    } else
        return std::move(worker.assignedStru);
}
