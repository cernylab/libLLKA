/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <QApplication>

#include "ui/main_window.h"

auto main(int argc, char **argv) -> int
{
    QApplication app{argc, argv};

    MainWindow mWin{};
    mWin.show();

    return app.exec();
}
