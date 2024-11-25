/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NTC_ASSIGNMENT_TABLE_MODEL
#define _NTC_ASSIGNMENT_TABLE_MODEL

#include <QAbstractTableModel>
#include <QItemDelegate>

#include "../../assignment.h"

class NtCAssignmentTableModel : public QAbstractTableModel {
    Q_OBJECT

public:
    explicit NtCAssignmentTableModel(QObject *parent = nullptr);

    auto data(const QModelIndex &index, int role = Qt::DisplayRole) const -> QVariant override;
    auto headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const -> QVariant override;

    auto columnCount(const QModelIndex &parent = {}) const -> int override;
    auto rowCount(const QModelIndex &parent = {}) const -> int override;

    auto structure() const -> const LoadedStructure &;
    auto setStructure(LoadedStructure structure) noexcept -> void;

private:
    LoadedStructure m_structure;
};

class NtCClassColorizerDelegate : public QItemDelegate {
public:
    explicit NtCClassColorizerDelegate(NtCAssignmentTableModel *model, QObject *parent = nullptr);

    auto paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const -> void override;

private:
    NtCAssignmentTableModel *h_model;
};

#endif // _NTC_ASSIGNMENT_TABLE_MODEL
