/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "main_window.h"
#include <ui_main_window.h>

#include "../assignment.h"
#include "../worker.h"
#include "assignment_parameters_dlg.h"
#include "detailed_step_info_dlg.h"
#include "setup_tracepoints_dlg.h"
#include "view_last_trace_dlg.h"
#include "models/ntc_assignment_table_model.h"

#include <llka_util.h>

#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QPushButton>
#include <QTableView>

#include <fstream>

inline const std::filesystem::path AVG_NU_ANGLES_FILE{"nu_angles.csv"};
inline const std::filesystem::path CLUSTERS_FILE{"clusters.csv"};
inline const std::filesystem::path CONFALS_FILE{"confals.csv"};
inline const std::filesystem::path GOLDEN_STEPS_FILE{"golden_steps.csv"};
inline const std::filesystem::path CONFAL_PERCENTILES_FILE{"confal_percentiles.csv"};

inline const LLKA_ClassificationLimits CLASSIFICATION_LIMITS{
    7,            // Minimum nearest neighbors
    11,           // Number of used nearest neighbors
    LLKA_deg2rad(28.0),    // Average neighbors torsion cutoff
    LLKA_deg2rad(28.0),    // Nearest neighbor torsion cutoff
    LLKA_deg2rad(60.0),    // Manhattan distance cutoff
    LLKA_deg2rad(72.0),    // Pseudorotation distance cutoff
    0.001111      // Minimum cluster votes
};

static
auto getClassificationContext(const std::filesystem::path &path, const AssignmentParameters &params) -> LLKA::ClassificationContext
{
    static const auto errMsg = QString{
        "Failed to load classification data for %5. "
        "Make sure that you set the correct path to the parameters directory and that the directory contains files %1, %2, %3 and %4"
    }.arg(AVG_NU_ANGLES_FILE.c_str(), CLUSTERS_FILE.c_str(), CONFALS_FILE.c_str(), GOLDEN_STEPS_FILE.c_str());

    auto resAvgNuAngles = LLKA::loadClusterNuAngles(path / AVG_NU_ANGLES_FILE);
    if (!resAvgNuAngles.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            errMsg.arg("average Nu angles"),
            QMessageBox::Ok
        );
        return {};
    }

    auto resClusters = LLKA::loadClusters(path / CLUSTERS_FILE);
    if (!resClusters.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            errMsg.arg("clusters"),
            QMessageBox::Ok
        );
        return {};
    }

    auto resConfals = LLKA::loadConfals(path / CONFALS_FILE);
    if (!resConfals.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            errMsg.arg("confals"),
            QMessageBox::Ok
        );
        return {};
    }

    auto resGoldenSteps = LLKA::loadGoldenSteps(path / GOLDEN_STEPS_FILE);
    if (!resGoldenSteps.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            errMsg.arg("golden steps"),
            QMessageBox::Ok
        );
        return {};
    }

    auto resConfalPercentiles = LLKA::loadConfalPercentiles(path / CONFAL_PERCENTILES_FILE);
    if (!resConfalPercentiles.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            errMsg.arg("confal percentiles"),
            QMessageBox::Ok
        );
        return {};
    }

    auto ctxRes = LLKA::initializeClassificationContext(
        resClusters.success(),
        resGoldenSteps.success(),
        resConfals.success(),
        resAvgNuAngles.success(),
        resConfalPercentiles.success(),
        CLASSIFICATION_LIMITS,
        params.maxCloseEnoughRmsd
    );
    if (!ctxRes.isSuccess()) {
        QMessageBox::warning(
            nullptr,
            "Cannot initialize classification context",
            QString{"Classification context initialization failed with error code %1"}.arg(LLKA::errorToString(ctxRes.failure()).c_str()),
            QMessageBox::Ok
        );
        return {};
    }

    return ctxRes.success();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow{parent},
    ui{new Ui::MainWindow}
{
    ui->setupUi(this);

    m_browseParams = new QFileDialog{this, "Select directory with classification parameters files", m_currentParamsPath, "CSV file (*.csv)"};
    m_browseParams->setFileMode(QFileDialog::Directory);

    m_browseCif = new QFileDialog{this, "Select mmCIF file with your structure", m_currentCifPath, "mmCIF file (*.cif *.mmcif)"};
    m_browseCif->setFileMode(QFileDialog::ExistingFile);

    m_setupTracepointsDlg = new SetupTracepointsDlg{this};

    m_assignmentModel = new NtCAssignmentTableModel{this};
    ui->qtbv_assignment->setModel(m_assignmentModel);

    auto ntcClassDelegate = new NtCClassColorizerDelegate{m_assignmentModel, this};
    ui->qtbv_assignment->setItemDelegateForColumn(2, ntcClassDelegate);

    m_params = AssignmentParameters::getDefault();

    connect(ui->qa_actionExit, &QAction::triggered, [this]() { close(); });
    connect(ui->qa_setupTracepoints, &QAction::triggered, [this]() {
        m_setupTracepointsDlg->readTracepointsState();
        m_setupTracepointsDlg->exec();
    });
    connect(ui->qa_viewLastTrace, &QAction::triggered, [this]() {
        ViewLastTraceDlg dlg{this};
        dlg.exec();
    });
    connect(ui->qa_assignmentParameters, &QAction::triggered, [this]() {
        AssignmentParametersDlg dlg{m_params, this};
        auto ret = dlg.exec();
        if (ret == QDialog::Accepted)
            this->m_params = dlg.params();
    });
    connect(ui->qtbv_assignment, &QTableView::doubleClicked, [this](const QModelIndex &index) {
        if (!index.isValid())
            return;

        auto row = index.row();
        auto col = index.column();

        if (row < 0 || row >= m_assignmentModel->rowCount() || col < 0 || col >= m_assignmentModel->columnCount())
            return;

        const auto &steps = m_assignmentModel->structure().steps;
        DetailedStepInfoDlg detail{steps, row, this};
        detail.exec();
    });
    connect(ui->qpb_browseParams, &QPushButton::clicked, [this]() {
        if (m_currentParamsPath.isEmpty())
            m_browseParams->setDirectory(QDir::currentPath());
        else
            m_browseParams->setDirectory(m_currentParamsPath);

        auto ret = m_browseParams->exec();
        if (ret == QDialog::Accepted) {
            m_currentParamsPath = m_browseParams->selectedFiles().constFirst();
            ui->qle_paramsDirPath->setText(m_currentParamsPath);
        }
    });
    connect(ui->qpb_browseCif, &QPushButton::clicked, [this]() {
        if (m_currentCifPath.isEmpty())
            m_browseCif->setDirectory(QDir::currentPath());
        else
            m_browseCif->setDirectory(QFileInfo{m_currentCifPath}.dir());

        auto ret = m_browseCif->exec();
        if (ret == QDialog::Accepted) {
            m_currentCifPath = m_browseCif->selectedFiles().constFirst();
            ui->qle_cifFilePath->setText(m_currentCifPath);
        }
    });
    connect(ui->qpb_assign, &QPushButton::clicked, this, &MainWindow::assign);
}

MainWindow::~MainWindow()
{
    delete ui;
}

auto MainWindow::assign() -> void
{
    if (m_currentParamsPath.isEmpty()) {
        QMessageBox::warning(
            this,
            "Empty path",
            "Path to the directory with classification parameters is empty.",
            QMessageBox::Ok
        );
        return;
    }
    std::string strParams = m_currentParamsPath.toUtf8().data();
    std::filesystem::path pathParams{strParams};

    if (m_currentCifPath.isEmpty()) {
        QMessageBox::warning(
            this,
            "Empty path",
            "Path to the structure file is empty.",
            QMessageBox::Ok
        );
        return;
    }
    std::string strCif = m_currentCifPath.toUtf8().data();
    std::filesystem::path pathCif{strCif};

    auto ctx = getClassificationContext(pathParams, m_params);
    if (!ctx.isValid())
        return;

    LLKA::trace(); // Clear previous trace

    auto assignedStructure = runAssignment(ctx, std::move(pathCif), this);
    if (assignedStructure.isValid()) {
        m_assignmentModel->setStructure(std::move(assignedStructure));
        ui->qtbv_assignment->resizeColumnsToContents();
    }

    auto traceOutputFile = m_setupTracepointsDlg->outputFile();
    if (!traceOutputFile.isEmpty()) {
        std::filesystem::path path{traceOutputFile.toUtf8().data()};
        std::ofstream out{path};

        out << LLKA::trace(true);
    }
}
