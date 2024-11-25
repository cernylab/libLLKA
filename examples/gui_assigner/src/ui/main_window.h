/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _MAIN_WINDOW_H
#define _MAIN_WINDOW_H

#include "../assignment_parameters.h"

#include <QMainWindow>

#include <vector>

class NtCAssignmentTableModel;
class QFileDialog;

class SetupTracepointsDlg;
class Step;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    auto assign() -> void;

    QFileDialog *m_browseParams;
    QString m_currentParamsPath;

    QFileDialog *m_browseCif;
    QString m_currentCifPath;

    NtCAssignmentTableModel *m_assignmentModel;
    SetupTracepointsDlg *m_setupTracepointsDlg;

    AssignmentParameters m_params;

    Ui::MainWindow *ui;
};

#endif // _MAIN_WINDOW_H
