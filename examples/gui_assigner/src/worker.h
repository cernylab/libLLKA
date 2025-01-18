/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _WORKER_H
#define _WORKER_H

#include "assignment.h"
#include "ui/assignment_in_progress_dlg.h"

#include <QObject>
#include <QString>

#include <atomic>
#include <filesystem>

class WorkerResult {
public:
    enum Outcome {
        Success,
        Failure,
        Canceled
    };

    Outcome outcome;
    QString message;

    WorkerResult(Outcome outcome, QString message) :
        outcome{outcome},
        message{std::move(message)}
    {
    }
};

class Worker : public QObject {
    Q_OBJECT

public:
    Worker(const LLKA::ClassificationContext &ctx, std::filesystem::path pathToCif);

    auto suceeded() const { return !m_cancel && !m_failed; }

    LoadedStructure assignedStru;

private:
    std::filesystem::path m_pathToCif;
    std::atomic<bool> m_cancel;
    bool m_failed;
    const LLKA::ClassificationContext &m_ctx;

public slots:
    void cancel();
    void work();

signals:
    void done(const WorkerResult &result);
    void reportProgress(AssignmentInProgressDlg::Stage stage, int percentDone);
};

auto runAssignment(const LLKA::ClassificationContext &ctx, std::filesystem::path pathToCif, QWidget *dialogParent) -> LoadedStructure;

#endif // _WORKER_H
