#pragma once

#include <QWidget>
#include "cgnslib.h"

class QTreeView;
class QStandardItem;
class BBase;
class Base;
class Zone;
class CgnsPanel;

class CgnsNode
{
public:
    CgnsNode(){};
    ~CgnsNode(){};
public:
    double id;
    double pid;
    std::string name;
    std::string label;
    std::string parent_name;
    std::string data_type;
    QStandardItem * item;
    std::string dimstr;
    cgsize_t np;
    cglong_t data_size;
    std::string valuestr;
    std::string link_file;
    std::string link_node;
};


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
    void ReadCgnsToQStandardItem(const QString &fileName);
    void ReadCgnsChild(int cgio_num, double pid, QStandardItem *pItem );
private:
    void onTreeItemClicked(const QModelIndex &index);
private:
    void AddCgnsToQStandardItem(BBase *bbase);
    void AddZoneItem(Base * base, QStandardItem *baseItem);
    void AddZoneSections(Zone * zone, QStandardItem *zoneItem);
    void AddZoneBcs(Zone * zone, QStandardItem *zoneItem);
    void AddZoneFlowSolution(Zone * zone, QStandardItem *zoneItem);
public:
    void SetCgnsPanel(CgnsPanel * cgnsPanel);
signals:
private:
    CgnsPanel * cgnsPanel = nullptr;
    QTreeView * treeView = nullptr;
    int fileId = -1;
    int baseId = -1;
};
