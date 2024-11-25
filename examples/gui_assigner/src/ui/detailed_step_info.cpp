/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "detailed_step_info.h"
#include <ui_detailed_step_info.h>

#include <llka_util.h>

#include <QStandardItemModel>

inline const QStringList METRICS_HEADER{
    LLKA::dinucleotideTorsionName(LLKA_TOR_DELTA_1, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_EPSILON_1, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_ZETA_1, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_ALPHA_2, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_BETA_2, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_GAMMA_2, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_DELTA_2, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_CHI_1, true).c_str(),
    LLKA::dinucleotideTorsionName(LLKA_TOR_CHI_2, true).c_str(),
    LLKA::crossResidueMetricName(LLKA_XR_DIST_CC, true).c_str(),
    LLKA::crossResidueMetricName(LLKA_XR_DIST_NN, true).c_str(),
    LLKA::crossResidueMetricName(LLKA_XR_TOR_MU, true).c_str()
};

inline const QStringList NU_ANGLES_HEADER{"ν_0", "ν_1", "ν_2", "ν_3", "ν_4"};

static
auto clearModel(QAbstractItemModel *model, bool columns)
{
    auto simodel = qobject_cast<QStandardItemModel *>(model);
    if (simodel) {
        if (columns)
            simodel->removeColumns(0, simodel->columnCount());
        else
            simodel->removeRows(0, simodel->rowCount());
    }
}

DetailedStepInfo::DetailedStepInfo(QWidget *parent) :
    QWidget{parent},
    ui{new Ui::DetailedStepInfo{}}
{
    ui->setupUi(this);
}

DetailedStepInfo::DetailedStepInfo(const Step &step, QWidget *parent) :
    QWidget{parent},
    ui{new Ui::DetailedStepInfo{}}
{
    ui->setupUi(this);

    setDetails(step);
}

DetailedStepInfo::~DetailedStepInfo()
{
    delete ui;
}

auto DetailedStepInfo::clearDetails() -> void
{
    ui->qle_assignedNtC->setText("");
    ui->qle_assignedCANA->setText("");
    ui->qle_closestNtC->setText("");
    ui->qle_closestCANA->setText("");
    ui->qle_rmsdClosest->setText("");
    ui->qle_closestGoldenStep->setText("");
    ui->qle_firstPseudorotation->setText("");
    ui->qle_secondPseudorotation->setText("");
    ui->qle_firstTau->setText("");
    ui->qle_secondTau->setText("");
    ui->ql_violationsList->setText("");
    ui->qle_firstSugarPucker->setText("");
    ui->qle_secondSugarPucker->setText("");

    clearModel(ui->qtbv_confalScores->model(), false);
    clearModel(ui->qtbv_firstNuAngles->model(), true);
    clearModel(ui->qtbv_secondNuAngles->model(), true);
    clearModel(ui->qtbv_firstNuAngleDifferences->model(), true);
    clearModel(ui->qtbv_secondNuAngleDifferences->model(), true);
}

auto DetailedStepInfo::fillNuAngles(QTableView *view, const LLKA_NuAngles &angles, bool asFullAngle) -> void
{
    QStandardItemModel *model = qobject_cast<QStandardItemModel *>(view->model());

    if (model == nullptr) {
        model = new QStandardItemModel{this};
        model->setVerticalHeaderLabels(NU_ANGLES_HEADER);

        view->setModel(model);
    } else
        model->removeColumns(0, model->columnCount());


    auto nu_0 = LLKA_rad2deg(angles.nu_0);
    auto nu_1 = LLKA_rad2deg(angles.nu_1);
    auto nu_2 = LLKA_rad2deg(angles.nu_2);
    auto nu_3 = LLKA_rad2deg(angles.nu_3);
    auto nu_4 = LLKA_rad2deg(angles.nu_4);

    if (asFullAngle) {
        nu_0 = LLKA_fullAngleFromDeg(nu_0);
        nu_1 = LLKA_fullAngleFromDeg(nu_1);
        nu_2 = LLKA_fullAngleFromDeg(nu_2);
        nu_3 = LLKA_fullAngleFromDeg(nu_3);
        nu_4 = LLKA_fullAngleFromDeg(nu_4);
    }

    model->setItem(0, 0, new QStandardItem{QString::number(nu_0)});
    model->setItem(1, 0, new QStandardItem{QString::number(nu_1)});
    model->setItem(2, 0, new QStandardItem{QString::number(nu_2)});
    model->setItem(3, 0, new QStandardItem{QString::number(nu_3)});
    model->setItem(4, 0, new QStandardItem{QString::number(nu_4)});

    view->resizeColumnsToContents();
}

auto DetailedStepInfo::setDetails(const Step &step) -> void
{
    switch (step.state) {
    case Step::Classified:
        ui->ql_stepNotAssignedReason->setVisible(false);
        break;
    case Step::NotYetClassified:
        ui->ql_stepNotAssignedReason->setText("This step has not been classified yet.");
        [[ fallthrough ]];
    case Step::ClassificationFailed:
        ui->ql_stepNotAssignedReason->setText("Classification of this step has failed.");
        ui->ql_stepNotAssignedReason->setVisible(true);
        clearDetails();
        return;
    }

    const auto &clsf = step.classification;

    ui->qle_assignedNtC->setText(LLKA::NtCToName(clsf.assignedNtC).c_str());
    ui->qle_assignedCANA->setText(LLKA::CANAToName(clsf.assignedCANA).c_str());

    ui->qle_closestNtC->setText(LLKA::NtCToName(clsf.closestNtC).c_str());
    ui->qle_closestCANA->setText(LLKA::CANAToName(clsf.closestCANA).c_str());

    ui->qle_rmsdClosest->setText(QString::number(clsf.rmsdToClosestNtC));

    ui->qle_closestGoldenStep->setText(clsf.closestGoldenStep.c_str());

    ui->qle_firstPseudorotation->setText(QString::number(LLKA_rad2deg(LLKA_fullAngleFromRad(clsf.ribosePseudorotation_1))));
    ui->qle_secondPseudorotation->setText(QString::number(LLKA_rad2deg(LLKA_fullAngleFromRad(clsf.ribosePseudorotation_2))));
    ui->qle_firstTau->setText(QString::number(LLKA_rad2deg(LLKA_fullAngleFromRad(clsf.tau_1))));
    ui->qle_secondTau->setText(QString::number(LLKA_rad2deg(LLKA_fullAngleFromRad(clsf.tau_2))));

    auto viols = clsf.namedViolations();
    QString violsStr{};
    for (const auto &v : viols)
        violsStr += QString{"- %1\n"}.arg(v.c_str());
    ui->ql_violationsList->setText(violsStr);

    ui->qle_firstSugarPucker->setText(LLKA::sugarPuckerToName(clsf.sugarPucker_1, LLKA_SPN_FANCY).c_str());
    ui->qle_secondSugarPucker->setText(LLKA::sugarPuckerToName(clsf.sugarPucker_2, LLKA_SPN_FANCY).c_str());

    if (!clsf.hasViolations()) {
        QStandardItemModel *model = qobject_cast<QStandardItemModel *>(ui->qtbv_confalScores->model());

        if (model == nullptr) {
            model = new QStandardItemModel{this};
            model->setHorizontalHeaderLabels(METRICS_HEADER);
            model->setHorizontalHeaderItem(12, new QStandardItem{"Total"});

            ui->qtbv_confalScores->setModel(model);
        } else
            model->removeRows(0, model->rowCount());

        model->setItem(0, 0, new QStandardItem{QString::number(clsf.confalScore.delta_1)});
        model->setItem(0, 1, new QStandardItem{QString::number(clsf.confalScore.epsilon_1)});
        model->setItem(0, 2, new QStandardItem{QString::number(clsf.confalScore.zeta_1)});
        model->setItem(0, 3, new QStandardItem{QString::number(clsf.confalScore.alpha_2)});
        model->setItem(0, 4, new QStandardItem{QString::number(clsf.confalScore.beta_2)});
        model->setItem(0, 5, new QStandardItem{QString::number(clsf.confalScore.gamma_2)});
        model->setItem(0, 6, new QStandardItem{QString::number(clsf.confalScore.delta_2)});
        model->setItem(0, 7, new QStandardItem{QString::number(clsf.confalScore.chi_1)});
        model->setItem(0, 8, new QStandardItem{QString::number(clsf.confalScore.chi_2)});
        model->setItem(0, 9, new QStandardItem{QString::number(clsf.confalScore.CC)});
        model->setItem(0, 10, new QStandardItem{QString::number(clsf.confalScore.NN)});
        model->setItem(0, 11, new QStandardItem{QString::number(clsf.confalScore.mu)});
        model->setItem(0, 12, new QStandardItem{QString::number(clsf.confalScore.total)});

        ui->qtbv_confalScores->resizeColumnsToContents();
    } else {
        auto model = qobject_cast<QStandardItemModel *>(ui->qtbv_confalScores->model());
        if (model)
            model->removeRows(0, model->rowCount());
    }

    fillNuAngles(ui->qtbv_firstNuAngles, clsf.nuAngles_1, true);
    fillNuAngles(ui->qtbv_secondNuAngles, clsf.nuAngles_2, true);
    fillNuAngles(ui->qtbv_firstNuAngleDifferences, clsf.nuAngleDifferences_1, false);
    fillNuAngles(ui->qtbv_secondNuAngleDifferences, clsf.nuAngleDifferences_2, false);
}
