/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "ntc_assignment_table_model.h"

#include "ntc_colors.h"

#include <QColor>
#include <QFontMetrics>
#include <QPainter>

#include <cmath>

/* https://alienryderflex.com/hsp.html */
auto luminance(const QColor &clr)
{
    auto r = clr.redF();
    auto g = clr.greenF();
    auto b = clr.blueF();

    return std::sqrt(0.299 * r * r + 0.587 * g * g + 0.114 * b * b);
}

static
auto rmsdToSemaphore(float rmsd) {
    float MinRmsd = 0.0f;
    float MaxRmsd = 1.0f;
    float Half = MaxRmsd / 2.0f;

    float normalized = (rmsd + MinRmsd) / MaxRmsd;
    if (normalized > MaxRmsd)
        normalized = MaxRmsd;
    int r = int((255.0f * (2.0f * normalized < 1.0f ? 2.0f * normalized : 1.0f)) + 0.5f);
    int g = int(255.0f * (1.0f - 2.0f * (normalized - Half > 0.0f ? normalized - Half : 0.0f)));

    return QColor(r, g, 0);
}

NtCClassColorizerDelegate::NtCClassColorizerDelegate(NtCAssignmentTableModel *model, QObject *parent) :
    QItemDelegate{parent},
    h_model{model}
{
}

auto NtCClassColorizerDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const -> void
{
    auto v = h_model->data(h_model->index(index.row(), 2), Qt::UserRole + 1);
    const auto ntc = qvariant_cast<LLKA_NtC>(v);
    const auto &clrs = NTC_COLORS.at(ntc);

    const auto ntcName = QString::fromStdString(LLKA::NtCToName(ntc));

    const auto w = option.rect.width();
    const auto h = option.rect.height();
    const auto fromX = option.rect.x();
    const auto fromY = option.rect.y();
    const auto wHalf = w / 2;

    painter->save();

    // First NtC color
    QBrush br{clrs.first};
    painter->fillRect(fromX, fromY, w, h, br);

    // Second NtC color
    br.setColor(clrs.second);
    painter->fillRect(fromX + wHalf, fromY, w - wHalf, h, br);

    // NtC name
    const auto textRect = option.fontMetrics.boundingRect(ntcName);
    const auto leftOffset = fromX + (w - textRect.width()) / 2;
    const auto topOffset = fromY + int(((h + textRect.height() / 2.0f) / 2.0f) + 0.5f);

    auto lum = luminance(clrs.first);
    painter->setPen(lum < 0.57 ? Qt::white : Qt::black);
    painter->drawText(leftOffset, topOffset, ntcName);

    painter->restore();
}

NtCAssignmentTableModel::NtCAssignmentTableModel(QObject *parent) :
    QAbstractTableModel{parent}
{
}

auto NtCAssignmentTableModel::columnCount(const QModelIndex &) const -> int
{
    return 5;
}

auto NtCAssignmentTableModel::data(const QModelIndex &index, int role) const -> QVariant
{
    auto row = index.row();
    auto col = index.column();

    if (row < 0 || row >= rowCount())
        return {};
    if (col < 0 || col >= columnCount())
        return {};

    const auto &step = m_structure.steps[row];

    if (role == Qt::BackgroundRole) {
        if (col == 4)
            return rmsdToSemaphore(step.classification.rmsdToClosestNtC);
    } else if (role == Qt::ForegroundRole) {
        if (col == 4) {
            auto clr = rmsdToSemaphore(step.classification.rmsdToClosestNtC);
            auto lum = luminance(clr);
            return lum < 0.57 ? QColor{Qt::white} : QColor{Qt::black};
        }
    } else if (role == Qt::DisplayRole) {
        switch (col) {
        case 0:
            return QString::fromStdString(step.name);
        case 1:
            return step.state == Step::Classified ? QString::fromStdString(LLKA::CANAToName(step.classification.assignedCANA)) : "-";
        case 2:
            return QString::fromStdString(LLKA::NtCToName(step.classification.assignedNtC));
        case 3:
            return step.state == Step::Classified ? QVariant{step.classification.confalScore.total} : QVariant{"-"};
        case 4:
            return step.state == Step::Classified ? QVariant{step.classification.rmsdToClosestNtC} : "-";
        }
    } else if (role == Qt::UserRole + 1) {
        switch (col) {
        case 2:
            return step.classification.assignedNtC;
        }
    }

    return {};
}

auto NtCAssignmentTableModel::headerData(int section, Qt::Orientation orientation, int role) const -> QVariant
{
    if (role != Qt::DisplayRole)
        return {};

    if (orientation == Qt::Vertical)
        return section + 1;
    else if (orientation == Qt::Horizontal) {
        switch (section) {
        case 0:
            return "Step name";
        case 1:
            return "CANA";
        case 2:
            return "NtC";
        case 3:
            return "Confal";
        case 4:
            return "RMSD";
        }

        return {};
    }

    return {};
}

auto NtCAssignmentTableModel::rowCount(const QModelIndex &) const -> int
{
    return m_structure.steps.size();
}

auto NtCAssignmentTableModel::setStructure(LoadedStructure structure) noexcept -> void
{
    beginResetModel();

    m_structure = std::move(structure);

    endResetModel();
}

auto NtCAssignmentTableModel::structure() const -> const LoadedStructure &
{
    return m_structure;
}
