/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "assignment_parameters_dlg.h"
#include <ui_assignment_parameters_dlg.h>

#include <QDialogButtonBox>

#include <cassert>

AssignmentParametersDlg::AssignmentParametersDlg(const AssignmentParameters &params, QWidget *parent) :
    QDialog{parent},
    ui{new Ui::AssignmentParametersDlg{}},
    m_params{params}
{
    ui->setupUi(this);

    assert(ui->qdspb_maxCloseEnoughRmsd->minimum() < m_params.maxCloseEnoughRmsd && m_params.maxCloseEnoughRmsd <= ui->qdspb_maxCloseEnoughRmsd->maximum());

    ui->qdspb_maxCloseEnoughRmsd->setValue(m_params.maxCloseEnoughRmsd);
    connect(ui->buttonBox, &QDialogButtonBox::accepted, [this]() {
        m_params.maxCloseEnoughRmsd = ui->qdspb_maxCloseEnoughRmsd->value();

        accept();
    });
    connect(ui->buttonBox, &QDialogButtonBox::rejected, [this]() { reject(); });
}

AssignmentParametersDlg::~AssignmentParametersDlg()
{
    delete ui;
}

auto AssignmentParametersDlg::params() const -> const AssignmentParameters &
{
    return m_params;
}
