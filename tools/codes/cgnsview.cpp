#include "cgnsview.h"
#include <QTreeView>
#include <QStandardItemModel>
#include <QHeaderView>
#include <QDir>
#include <QDirIterator>
#include <QResizeEvent>
#include "cgnslib.h"
#include "CgnsBase.h"
#include "CgnsBc.h"
#include <iostream>

CgnsView::CgnsView(QWidget *parent)
    : QWidget{parent}
{
    this->treeView = new QTreeView(this);
    this->treeView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

    this->SetUpCgnsInfo();
    QObject::connect(this->treeView, &QAbstractItemView::clicked, this, &CgnsView::onTreeItemClicked);
}

void CgnsView::resizeEvent(QResizeEvent *event)
{
    this->treeView->resize(event->size());

}

void CgnsView::SetUpHeader()
{
    QStandardItemModel *itemModel = new QStandardItemModel();
    this->treeView->setModel(itemModel);
    this->treeView->header()->setSectionResizeMode(QHeaderView::Stretch);
    this->treeView->show();
    QStringList labels;
    labels << "CGNSView";
    itemModel->setHorizontalHeaderLabels(labels);
}

void printItem(QStandardItem* item, int depth = 0)
{
    if (!item) return; // 如果 item 为空，返回

    // 打印当前项的文本
    qDebug().nospace() << QString("%1%2").arg(QString("  ").repeated(depth)).arg(item->text());

    // 遍历子项
    for (int i = 0; i < item->rowCount(); ++i) {
        printItem(item->child(i), depth + 1); // 递归调用
    }
}

void CgnsView::SetUpCgnsInfo()
{
    this->SetUpHeader();

    QStandardItemModel *itemModel = static_cast<QStandardItemModel *>(this->treeView->model());

    qDebug() << " dataTypeName 1= " << DataTypeName[ RealDouble ];
    qDebug() << " dataTypeName 2= " << DataTypeName[ ComplexSingle ];

    QString fileName = "D:/work/cgns_work/cgnsgrid_example/yf17_hdf5.cgns";

    // QStandardItem *rootItem = new QStandardItem( "/" );
    // itemModel->appendRow(rootItem);

    QModelIndex rootIndex = this->treeView->rootIndex();
    qDebug() << "rootIndex = " << rootIndex;

    int count = itemModel->rowCount();
    qDebug() << "rowCount=" << count;

    // 打印每个项的信息
    for (int i = 0; i < itemModel->rowCount(); ++i)
    {
        printItem(itemModel->item(i));
    }

    int row = rootIndex.row();
    int col = rootIndex.column();
    qDebug() <<"rootIndex row=" << row << "column=" << col;

    QString text = rootIndex.data(Qt::DisplayRole).toString();
    qDebug() << "rootIndex text = " << text;

    this->treeView->expandAll();

    ::ReadCgnsFile( fileName.toLatin1().data() );

    this->AddCgnsToQStandardItem( ::GetBBase() );
}

void CgnsView::AddCgnsToQStandardItem(BBase *bbase)
{
    this->fileId = bbase->fileId;
    qDebug() << "";
    qDebug() << "bbase->nBases = " << bbase->nBases;
    qDebug() << "bbase->fileId = " << bbase->fileId;
    qDebug() << "bbase->precision = " << bbase->precision;

    qDebug() << "   CGNS File Type = " << bbase->file_type << " FileTypeName = " << ::GetCgnsFileTypeName( bbase->file_type );

    QHeaderView* header = this->treeView->header(); // 获取 QHeaderView
    qDebug() << "Header Count:" << header->count(); // 打印列的数量

    // 打印每一列的标签
    for (int section = 0; section < header->count(); ++section)
    {
        qDebug() << "Header Section" << section << ":" << header->model()->headerData(section, Qt::Horizontal);
    }

    QModelIndex rootIndex = this->treeView->rootIndex();
    qDebug() << "rootIndex = " << rootIndex;
    int row = rootIndex.row();
    int col = rootIndex.column();
    qDebug() <<"rootIndex row=" << row << "column=" << col;

    QString text = rootIndex.data(Qt::DisplayRole).toString();
    qDebug() << "rootIndex text = " << text;

    QStandardItemModel *itemModel = static_cast<QStandardItemModel *>(this->treeView->model());
    //qDebug() <<"itemModel->item(0)->text() = " << itemModel->item(0)->text();
    //QStandardItem *rootItem = new QStandardItem( itemModel->item(0)->text() );

    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");

    QStandardItem *rootItem = new QStandardItem( "/" );
    itemModel->appendRow(rootItem);
    rootItem->setIcon(folderIcon);

    QStandardItem *libraryVersionItem = new QStandardItem("CGNSLibraryVersion");
    itemModel->appendRow(libraryVersionItem);
    libraryVersionItem->setIcon(fileIcon);

    this->treeView->expandRecursively(rootItem->index());

    for ( int iBase = 0; iBase < bbase->nBases; ++ iBase )
    {
        int baseId = iBase + 1;
        this->baseId = baseId;
        //std::vector<Base *> bases;
        Base * base = bbase->bases[ iBase ];
        qDebug() << "base->name=" << base->name;
        QStandardItem *baseItem = new QStandardItem( base->name );
        itemModel->appendRow(baseItem);
        baseItem->setIcon(folderIcon);

        this->AddZoneItem(base, baseItem);

        QStandardItem *dataClassItem = new QStandardItem( "DataClass" );
        baseItem->appendRow(dataClassItem);
        dataClassItem->setIcon(fileIcon);

        QStandardItem *dimensionalUnitsItem = new QStandardItem( "DimensionalUnits" );
        baseItem->appendRow(dimensionalUnitsItem);
        dimensionalUnitsItem->setIcon(fileIcon);

        this->treeView->expand(baseItem->index());
    }
}

void CgnsView::AddZoneItem(Base * base, QStandardItem *baseItem)
{
    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");
    qDebug() << "base->nZones=" << base->nZones;
    for ( int iZone = 0; iZone < base->nZones; ++ iZone )
    {
        qDebug() << "iZone=" << iZone;
        Zone * zone = base->zones[iZone];
        qDebug() << "zone->zoneName=" << zone->zoneName;
        QStandardItem *zoneItem = new QStandardItem( zone->zoneName );
        baseItem->appendRow(zoneItem);
        zoneItem->setIcon(folderIcon);

        QStandardItem *zoneTypeItem = new QStandardItem( "ZoneType" );
        zoneItem->appendRow(zoneTypeItem);
        zoneTypeItem->setIcon(fileIcon);

        qDebug() << "zone->zoneType=" << zone->zoneType;

        QStandardItem *zoneCoorItem = new QStandardItem( "GridCoordinates" );
        zoneItem->appendRow(zoneCoorItem);
        zoneCoorItem->setIcon(folderIcon);

        this->treeView->expandRecursively(zoneItem->index());

        for ( int iCoord = 0; iCoord < zone->nCoords; ++ iCoord )
        {
            Coor * coor = zone->coors[iCoord];
            qDebug() << "coor->coorName=" << coor->coorName;
            QStandardItem *coorItem = new QStandardItem( coor->coorName );
            zoneCoorItem->appendRow(coorItem);
            coorItem->setIcon(fileIcon);
        }

        this->AddZoneSections(zone, zoneItem);
        this->AddZoneFlowSolution(zone, zoneItem);
        this->AddZoneBcs(zone, zoneItem);
    }
}

void CgnsView::AddZoneSections(Zone * zone, QStandardItem *zoneItem)
{
    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");
    qDebug() << "zone->nSections=" << zone->nSections;
    for ( int iSection = 0; iSection < zone->nSections; ++ iSection )
    {
        Section * section = zone->sections[iSection];
        qDebug() << "section->sectionName=" << section->sectionName;
        QStandardItem *sectionItem = new QStandardItem( section->sectionName );
        zoneItem->appendRow(sectionItem);
        sectionItem->setIcon(folderIcon);

        QStandardItem *elemRangeItem = new QStandardItem( "ElementRange" );
        sectionItem->appendRow(elemRangeItem);
        elemRangeItem->setIcon(fileIcon);

        QStandardItem *elemConnItem = new QStandardItem( "ElementConnectivity" );
        sectionItem->appendRow(elemConnItem);
        elemConnItem->setIcon(fileIcon);

        QStandardItem *elemParentDataItem = new QStandardItem( "ParentData" );
        sectionItem->appendRow(elemParentDataItem);
        elemParentDataItem->setIcon(fileIcon);
    }
}

void CgnsView::AddZoneFlowSolution(Zone * zone, QStandardItem *zoneItem )
{
    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");

    std::cout << " CgnsView::AddZoneFlowSolution " << std::endl;

    std::cout << " this->fileId = " << this->fileId << std::endl;
    std::cout << " this->baseId = " << this->baseId << std::endl;
    std::cout << " zone->zoneId = " << zone->zoneId << std::endl;

    for ( int iSolution = 0; iSolution < zone->nSolutions; ++ iSolution )
    {
        int solutionId = iSolution + 1;
        Solution * solution = zone->solutions[iSolution];
        std::cout << " solutionName = " << solution->solutionName << std::endl;

        QStandardItem *solutionItem = new QStandardItem( solution->solutionName );
        zoneItem->appendRow(solutionItem);
        solutionItem->setIcon(folderIcon);

        for ( int iField = 0; iField < solution->nFields; ++ iField )
        {
            int fieldId = iField + 1;
            Field * field = solution->fields[iField];

            QStandardItem *fieldItem = new QStandardItem( field->fieldName );
            solutionItem->appendRow(fieldItem);
            if ( field->dimflag < 0 )
            {
                fieldItem->setIcon(fileIcon);
            }
            else
            {
                fieldItem->setIcon(folderIcon);
                if( field->dimflag == 0 )
                {
                    QStandardItem *classItem = new QStandardItem( "DataClass" );
                    fieldItem->appendRow(classItem);
                    classItem->setIcon(fileIcon);
                }
                else
                {
                    QStandardItem *dimItem = new QStandardItem( "DimensionalExponents" );
                    fieldItem->appendRow(dimItem);
                    dimItem->setIcon(fileIcon);
                }
            }
        }
   }

}

void CgnsView::AddZoneBcs(Zone * zone, QStandardItem *zoneItem)
{
    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");
    Bc * bc = zone->bc;
    qDebug() << "bc->nBocos=" << bc->nBocos;

    QStandardItem *zoneBcItem = new QStandardItem( "ZoneBC" );
    zoneItem->appendRow(zoneBcItem);
    zoneBcItem->setIcon(folderIcon);

    for ( int iBoco = 0; iBoco < bc->nBocos; ++ iBoco )
    {
        qDebug() << "\n";
        qDebug() <<  "-->iBoco  = " << iBoco << " bc->nBocos = " << bc->nBocos;
        int bocoId = iBoco + 1;
        Boco * boco = bc->bocos[ iBoco ];
        qDebug() << "boco->name=" << boco->name;
        QStandardItem *bcItem = new QStandardItem( boco->name.c_str() );
        zoneBcItem->appendRow(bcItem);
        bcItem->setIcon(folderIcon);

        QStandardItem *pointRangeItem = new QStandardItem( "PointRange" );
        bcItem->appendRow(pointRangeItem);
        pointRangeItem->setIcon(fileIcon);

        QStandardItem *gridLocationItem = new QStandardItem( "GridLocation" );
        bcItem->appendRow(gridLocationItem);
        gridLocationItem->setIcon(fileIcon);
    }
}

void CgnsView::onTreeItemClicked(const QModelIndex &index)
{
    qDebug() <<"CgnsView::on_treeView_clicked" << "index=" << index;
    int row = index.row();
    int col = index.column();
    qDebug() <<"item row=" << row << "column=" << col;

    QString text = index.data(Qt::DisplayRole).toString();
    qDebug() << "text = " << text;

}

