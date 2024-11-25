/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _ASSIGNMENT_PROGRESS_DLG_H
#define _ASSIGNMENT_PROGRESS_DLG_H

#include <QDialog>

namespace Ui {
    class AssignmentInProgressDlg;
}

class AssignmentInProgressDlg : public QDialog {
    Q_OBJECT

public:
    enum Stage {
        ReadingCif,
        Assigning
    };

    explicit AssignmentInProgressDlg(QWidget *parent = nullptr);
    ~AssignmentInProgressDlg();

private:
    Ui::AssignmentInProgressDlg *ui;

public slots:
    void onDone();
    void onProgress(Stage stage, int percentDone);

signals:
    void cancel();
};

#endif // _ASSIGNMENT_PROGRESS_DLG_H
