/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "detailed_step_info_dlg.h"
#include <ui_detailed_step_info_dlg.h>

#include "detailed_step_info.h"

QSize DetailedStepInfoDlg::m_lastDlgSize{};

DetailedStepInfoDlg::DetailedStepInfoDlg(const std::vector<Step> &steps, int stepIdx, QWidget *parent) :
    QDialog{parent},
    ui{new Ui::DetailedStepInfoDlg{}},
    m_steps{steps}
{
    ui->setupUi(this);

    for (size_t idx = 0; idx < steps.size(); idx++) {
        const auto &step = steps[idx];
        ui->qcbox_steps->addItem(QString::fromStdString(step.name), QVariant{qulonglong(idx)});
    }

    connect(ui->qcbox_steps, &QComboBox::activated, [this](int idx){
        auto stepIdx = ui->qcbox_steps->itemData(idx);
        if (m_details)
            m_details->setDetails(m_steps[stepIdx.toULongLong()]);
    });

    if (steps.size() > size_t(stepIdx)) {
        m_details = new DetailedStepInfo{steps[stepIdx], this};
        ui->qlay_main->addWidget(m_details);
        ui->qcbox_steps->setCurrentIndex(stepIdx);
    }
}

DetailedStepInfoDlg::~DetailedStepInfoDlg()
{
    delete ui;
}

auto DetailedStepInfoDlg::exec() -> int
{
    if (m_lastDlgSize.isValid())
        resize(m_lastDlgSize);

    auto ret = QDialog::exec();
    m_lastDlgSize = size();

    return ret;
}
