/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _SETUP_TRACEPOINTS_DLG_H
#define _SETUP_TRACEPOINTS_DLG_H

#include <QDialog>

#include <map>

class QCheckBox;
class QFileDialog;

namespace Ui {
    class SetupTracepointsDlg;
}

class SetupTracepointsDlg : public QDialog {
    Q_OBJECT

public:
    SetupTracepointsDlg(QWidget *parent = nullptr);
    ~SetupTracepointsDlg();

    auto outputFile() const -> QString;
    auto readTracepointsState() -> void;

private:
    Ui::SetupTracepointsDlg *ui;

    std::map<int32_t, bool> m_tracepointsToToggle;
    std::vector<QCheckBox *> m_checkboxes;

    QFileDialog *m_fileDlg;
};

#endif // _SETUP_TRACEPOINTS_DLG_H
