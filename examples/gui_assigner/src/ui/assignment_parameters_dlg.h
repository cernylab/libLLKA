/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _ASSIGNMENT_PARAMETERS_DLG_H
#define _ASSIGNMENT_PARAMETERS_DLG_H

#include "../assignment_parameters.h"

#include <QDialog>

namespace Ui {
    class AssignmentParametersDlg;
}

class AssignmentParametersDlg : public QDialog {
    Q_OBJECT

public:
    explicit AssignmentParametersDlg(const AssignmentParameters &params, QWidget *parent = nullptr);
    ~AssignmentParametersDlg();

    auto params() const -> const AssignmentParameters &;

private:
    Ui::AssignmentParametersDlg *ui;

    AssignmentParameters m_params;
};

#endif //  _ASSIGNMENT_PARAMETERS_DLG_H
