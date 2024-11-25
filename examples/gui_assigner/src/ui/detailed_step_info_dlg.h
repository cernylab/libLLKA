/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _DETAILED_STEP_INFO_DLG_H
#define _DETAILED_STEP_INFO_DLG_H

#include <QDialog>

#include "../assignment.h"

class DetailedStepInfo;

namespace Ui {
    class DetailedStepInfoDlg;
}

class DetailedStepInfoDlg : public QDialog {
    Q_OBJECT

public:
    DetailedStepInfoDlg(const std::vector<Step> &steps, int stepIdx, QWidget *parent = nullptr);
    ~DetailedStepInfoDlg();

    auto exec() -> int override;

private:
    Ui::DetailedStepInfoDlg *ui;

    DetailedStepInfo *m_details;
    const std::vector<Step> &m_steps;

    static QSize m_lastDlgSize;
};

#endif // _DETAILED_STEP_INFO_DLG_H

