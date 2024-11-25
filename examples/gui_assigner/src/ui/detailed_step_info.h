/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _DETAILED_STEP_INFO_H
#define _DETAILED_STEP_INFO_H

#include <QWidget>

#include "../assignment.h"

class QTableView;

namespace Ui {
    class DetailedStepInfo;
}

class DetailedStepInfo : public QWidget {
    Q_OBJECT

public:
    explicit DetailedStepInfo(QWidget *parent = nullptr);
    DetailedStepInfo(const Step &step, QWidget *parent);
    ~DetailedStepInfo();

    auto setDetails(const Step &step) -> void;

private:
    auto clearDetails() -> void;
    auto fillNuAngles(QTableView *view, const LLKA_NuAngles &angles, bool asFullAngle) -> void;

    Ui::DetailedStepInfo *ui;
};

#endif // _DETAILED_STEP_INFO_H
