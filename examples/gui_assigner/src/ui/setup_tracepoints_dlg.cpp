/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "setup_tracepoints_dlg.h"
#include <ui_setup_tracepoints_dlg.h>

#include <llka_cpp.h>

#include <QBoxLayout>
#include <QCheckBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QSpacerItem>

SetupTracepointsDlg::SetupTracepointsDlg(QWidget *parent) :
    QDialog{parent},
    ui{new Ui::SetupTracepointsDlg{}}
{
    ui->setupUi(this);

    m_fileDlg = new QFileDialog{nullptr};
    m_fileDlg->setAcceptMode(QFileDialog::AcceptSave);

    Q_ASSERT(ui->qscrl_main->widget() != nullptr);

    auto w = ui->qscrl_main->widget();
    auto lay = qobject_cast<QBoxLayout *>(w->layout());
    Q_ASSERT(lay != nullptr);

    ui->qle_filePath->setEnabled(false);

    const auto tpInfos = LLKA::tracepointInfo();
    m_checkboxes.reserve(tpInfos.size());
    for (const auto &tpi : tpInfos) {
        auto TPID = tpi.TPID;
        auto isEnabled = LLKA::tracepointState(TPID);
        auto cb = new QCheckBox{tpi.description.c_str(), this};
        cb->setCheckState(isEnabled ? Qt::Checked : Qt::Unchecked);
        connect(cb, &QCheckBox::checkStateChanged, [this, TPID](int state) {
            const bool enable = state == Qt::Checked;
            m_tracepointsToToggle[TPID] = enable;
        });

        cb->setProperty("TPID", TPID);
        m_checkboxes.push_back(cb);

        lay->addWidget(cb);
    }
    lay->addSpacerItem(new QSpacerItem{1, 1, QSizePolicy::Expanding, QSizePolicy::Expanding});

    connect(ui->qcb_outputFile, &QCheckBox::checkStateChanged, [this](int state) {
        const bool enable = state == Qt::Checked;

        ui->qle_filePath->setEnabled(enable);
        ui->qpb_browse->setEnabled(enable);
    });
    connect(ui->qpb_browse, &QPushButton::clicked, [this]() {
        if (!ui->qle_filePath->text().isEmpty())
            m_fileDlg->setDirectory(QFileInfo{ui->qle_filePath->text()}.absolutePath());

        auto ret = m_fileDlg->exec();
        if (ret == QDialog::Accepted)
            ui->qle_filePath->setText(m_fileDlg->selectedFiles().constFirst());
    });
    connect(ui->buttonBox, &QDialogButtonBox::accepted, [this]() {
        for (const auto &[TPID, state] : m_tracepointsToToggle)
            LLKA::toggleTracepoint(TPID, state);

        m_tracepointsToToggle.clear();
        accept();
    });
    connect(ui->buttonBox, &QDialogButtonBox::rejected, [this]() { reject(); });
}

SetupTracepointsDlg::~SetupTracepointsDlg()
{
    delete ui;
}

auto SetupTracepointsDlg::outputFile() const -> QString
{
    if (ui->qcb_outputFile->isChecked())
        return ui->qle_filePath->text();
    return {};
}

auto SetupTracepointsDlg::readTracepointsState() -> void
{
    for (auto &cb : m_checkboxes) {
        auto TPID = qvariant_cast<int32_t>(cb->property("TPID"));
        auto isEnabled = LLKA::tracepointState(TPID);
        cb->setCheckState(isEnabled ? Qt::Checked : Qt::Unchecked);
    }
}
