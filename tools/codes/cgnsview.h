#pragma once

#include <QWidget>
class QTreeView;
class QStandardItem;
class BBase;
class Base;
class Zone;

class CgnsView : public QWidget
{
    Q_OBJECT
public:
    explicit CgnsView(QWidget *parent = nullptr);
protected:
    void resizeEvent(QResizeEvent *event) override;
private:
    void SetUpCgnsInfo();
    void SetUpHeader();
private:
    void onTreeItemClicked(const QModelIndex &index);
private:
    void AddCgnsToQStandardItem(BBase *bbase);
    void AddZoneItem(Base * base, QStandardItem *baseItem);
    void AddZoneSections(Zone * zone, QStandardItem *zoneItem);
    void AddZoneBcs(Zone * zone, QStandardItem *zoneItem);
    void AddZoneFlowSolution(Zone * zone, QStandardItem *zoneItem);
signals:
private:
    QTreeView * treeView = nullptr;
    int fileId = -1;
    int baseId = -1;
};
