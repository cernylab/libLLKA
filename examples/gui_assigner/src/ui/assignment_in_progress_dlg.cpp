/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "assignment_in_progress_dlg.h"
#include <ui_assignment_in_progress_dlg.h>

#include <iostream>

AssignmentInProgressDlg::AssignmentInProgressDlg(QWidget *parent) :
    QDialog{parent},
    ui{new Ui::AssignmentInProgressDlg{}}
{
    ui->setupUi(this);

    connect(ui->buttonBox, &QDialogButtonBox::rejected, [this]() { emit cancel(); });
}

AssignmentInProgressDlg::~AssignmentInProgressDlg()
{
    delete ui;
}

void AssignmentInProgressDlg::onDone()
{
    accept();
}

void AssignmentInProgressDlg::onProgress(Stage stage, int percentDone)
{
    if (stage == ReadingCif)
        ui->ql_progress->setText("Reading mmCIF...");
    else if (stage == Assigning)
        ui->ql_progress->setText(QString{"Assigning: %1 % done..."}.arg(percentDone));
}
