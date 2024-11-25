/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "view_last_trace_dlg.h"

#include <llka_cpp.h>

#include <QDialogButtonBox>
#include <QTextEdit>
#include <QVBoxLayout>

ViewLastTraceDlg::ViewLastTraceDlg(QWidget *parent) :
    QDialog{parent}
{
    auto lay = new QVBoxLayout{this};
    auto text = new QTextEdit{this};
    text->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    text->setReadOnly(true);
    text->setAcceptRichText(false);

    auto bbox = new QDialogButtonBox{this};
    bbox->addButton(QDialogButtonBox::Cancel);

    lay->addWidget(text);
    lay->addWidget(bbox);

    auto trace = LLKA::trace(true);
    text->setPlainText(QString::fromStdString(trace));

    connect(bbox, &QDialogButtonBox::rejected, [this]() { accept(); });
}
