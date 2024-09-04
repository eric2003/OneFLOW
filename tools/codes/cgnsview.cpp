#include "cgnsview.h"
#include "cgnspanel.h"
#include <QTreeView>
#include <QStandardItemModel>
#include <QHeaderView>
#include <QDir>
#include <QDirIterator>
#include <QResizeEvent>
#include "cgnslib.h"
#include "CgnsBase.h"
#include "CgnsBc.h"
#include "cgns_io.h"
#include <iostream>
#include <format>
#include <map>

static int cgioNum = 0;
static double RootID = 0.0;
std::map<double,CgnsNode *> nodeMap;
std::map<QStandardItem *,CgnsNode *> cgnsNodeMap;

#define I4 int
#define U4 unsigned int
#define I8 cglong_t
#define U8 cgulong_t
#define R4 float
#define R8 double
#define X4 float
#define X8 double
#define C1 char
#define B1 unsigned char

enum DataTags {
    MTdata = 0,
    I4data,
    I8data,
    U4data,
    U8data,
    R4data,
    R8data,
    X4data,
    X8data,
    C1data,
    B1data,
    LKdata
};

static struct DataType {
    const char *name;
    int type;
    int bytes;
    int size;
} DataList[] = {
    {"MT", MTdata,  0, 0},
    {"I4", I4data,  4, sizeof(I4)},
    {"I8", I8data,  8, sizeof(I8)},
    {"U4", U4data,  4, sizeof(U4)},
    {"U8", U8data,  8, sizeof(U8)},
    {"R4", R4data,  4, sizeof(R4)},
    {"R8", R8data,  8, sizeof(R8)},
    {"X4", X4data,  8, sizeof(X4)*2},
    {"X8", X8data, 16, sizeof(X8)*2},
    {"C1", C1data,  1, 1},
    {"B1", B1data,  1, 1},
    {"LK", LKdata,  0, 0}
};

cgsize_t numberOfElements(int cgio_num, double nodeid)
{
    int n, ndim;
    cgsize_t dims[CGIO_MAX_DIMENSIONS];

    if ( cgio_get_dimensions (cgioNum, nodeid, &ndim, dims) ) {
        std::cout << "cgio_get_dimensions error\n";
        return -1;
    }

    if ( ndim == 0 )
    {
        return 0;
    }

    cgsize_t np = 1;

    for (int n = 0; n < ndim; n++)
    {
        np *= dims[n];
    }
    return np;
}

template<typename T>
void nestedFormatString( std::vector<char> &values, int ndim, cgsize_t * dims, int nElements, int depth, std::vector<int> &indices, int &icount, std::string & valuestr )
{
    if ( nElements > 200 ) return;

    if ( depth == ndim )
    {
        const T *x = (T *)values.data();

        double number = x[icount++];

        //valuestr += std::to_string(x[icount++]);
        valuestr += std::format("{:g}", number);
        valuestr += " ";
        return;
    }

    cgsize_t N = dims[ depth ];
    for ( int i = 0; i < N; ++i )
    {
        indices.push_back(i);
        nestedFormatString<T>( values, ndim, dims, nElements, depth + 1, indices, icount, valuestr );
        indices.pop_back();
        if( i == ( N - 1 ) && depth == ( ndim - 1 ) )
        {
            valuestr += "\n";
        }
    }
}

template<typename T>
void nestedFormatStringComplex( std::vector<char> &values, int ndim, cgsize_t * dims, int nElements, int depth, std::vector<int> &indices, int &icount, std::string & valuestr )
{
    if ( nElements > 200 ) return;

    if ( depth == ndim )
    {
        T *x = (T *)values.data();

        valuestr += "{";
        valuestr += std::to_string(x[icount++]);
        valuestr += ",";
        valuestr += std::to_string(x[icount++]);
        valuestr += "} ";

        return;
    }

    cgsize_t N = dims[ depth ];
    for ( int i = 0; i < N; ++i )
    {
        indices.push_back(i);
        nestedFormatStringComplex<T>( values, ndim, dims, nElements, depth + 1, indices, icount, valuestr );
        indices.pop_back();
        if( i == ( N - 1 ) && depth == ( ndim - 1 ) )
        {
            valuestr += "\n";
        }
    }
}


void printNestedFormatValues( const char *data_type, std::vector<char> &values, int ndim, cgsize_t * dims, int nElements, std::string & valuestr )
{
    int icount = 0;
    std::vector<int> indices = {};
    int depth = 0;

    qDebug() << "printValues data_type = " << data_type;
    if ( std::strncmp(data_type,"C1",2) == 0 ||
        std::strncmp(data_type,"B1",2) == 0 )
    {
        qDebug() << "std::strncmp(data_type,\"C1\",2) data_type = " << data_type;
        valuestr += values.data();
        qDebug() << "valuestr=" << valuestr;
    }
    else if ( std::strncmp(data_type,"I4",2) == 0 )
    {
        nestedFormatString<I4>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"I8",2) == 0 )
    {
        nestedFormatString<I8>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"U4",2) == 0 )
    {
        nestedFormatString<U4>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"U8",2) == 0 )
    {
        nestedFormatString<U8>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"R4",2) == 0 )
    {
        nestedFormatString<R4>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"R8",2) == 0 )
    {
        nestedFormatString<R8>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"X4",2) == 0 )
    {
        nestedFormatString<X4>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
    else if ( std::strncmp(data_type,"X8",2) == 0 )
    {
        nestedFormatString<X8>( values, ndim, dims, nElements, depth, indices, icount, valuestr );
    }
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

void CgnsView::SetUpCgnsInfo()
{
    this->SetUpHeader();

    QString fileName = "D:/work/cgns_work/cgnsgrid_example/yf17_hdf5.cgns";
    //QString fileName = "D:/work/cgns_work/ModernCGNS/codes/UserGuideCode/read_complex/complex_data.cgns";
    //QString fileName = "D:/work/cgns_work/ModernCGNS/codes/UserGuideCode/cg_link_write/grid/links.cgns";

    ::ReadCgnsFile( fileName.toLatin1().data() );

    //this->AddCgnsToQStandardItem( ::GetBBase() );
    this->ReadCgnsToQStandardItem( fileName );
}

void CgnsView::ReadCgnsChild(int cgio_num, double pid, QStandardItem *pItem )
{
    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");
    QIcon forwardIcon(":/forward.png");

    char parent_name[33];
    cgio_get_name(cgio_num, pid, parent_name);
    qDebug() << "parent_name = " << parent_name;

    char parent_label[33];
    cgio_get_label (cgio_num, pid, parent_label);
    qDebug() << "parent_label = " << parent_label;

    int num_children=-1;
    cgio_number_children (cgio_num, pid, &num_children);
    qDebug() << "num_children = " << num_children;
    if ( num_children <= 0 ) return;

    double childID = -1;
    char childName[33];
    for ( int n = 1; n <= num_children; n++ )
    {
        int cnt=-1;
        cgio_children_ids(cgio_num, pid, n, 1, &cnt, &childID);
        qDebug() << " n= " << n << " cnt = " << cnt;
        qDebug("childID = %30.22e", childID);
        cgio_get_name(cgio_num, childID, childName);
        qDebug() << "parent_name = " << parent_name;
        qDebug() << "childName = " << childName;
        char pid_name[33];
        cgio_get_name(cgio_num, pid, pid_name);
        qDebug() << "pid_name = " << pid_name;

        double node_id=-1;
        cgio_get_node_id (cgio_num, pid, childName, &node_id);
        qDebug("pid = %30.22e", pid);
        qDebug("node_id = %30.22e", node_id);
        qDebug("node_id = %f", node_id);
        qDebug("childID = %e", childID);

        int link_len = -1;
        cgio_is_link (cgio_num, node_id, &link_len);
        qDebug() << "link_len = " << link_len;

        std::string link_file;
        std::string link_node;

        if ( link_len > 0 )
        {
            int file_len = -1;
            int name_len = -1;
            cgio_link_size (cgio_num, node_id, &file_len, &name_len);
            qDebug() << "file_len = " << file_len << " name_len = " << name_len;
            std::vector<char> filename_vchar(file_len);
            std::vector<char> name_in_file_vchar(name_len);
            //cgio_get_link (cgio_num, node_id, filename, name_in_file);
            cgio_get_link (cgio_num, node_id, filename_vchar.data(), name_in_file_vchar.data());
            link_file = filename_vchar.data();
            link_node = name_in_file_vchar.data();

            qDebug() << "link_file = " << link_file << " link_node = " << link_node;
        }

        char label[CGIO_MAX_LABEL_LENGTH+1];
        cgio_get_label (cgio_num, node_id, label);
        qDebug() << "label = " << label;

        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type (cgio_num, node_id, data_type);
        qDebug() << "data_type = " << data_type;

        cglong_t data_size;
        cgio_get_data_size(cgio_num, node_id, &data_size);
        qDebug() << "data_size = " << data_size;

        int ndim;
        cgsize_t dims[CGIO_MAX_DIMENSIONS];
        cgio_get_dimensions (cgio_num, node_id, &ndim, dims);
        qDebug() << "ndim = " << ndim;
        std::string dimstr;

        qDebug() << "dims = ";
        for( int i = 0; i < ndim; ++ i )
        {
            qDebug() << dims[i];
            dimstr += std::to_string(dims[i]);
            if ( i != ndim-1 )
            {
                dimstr += " ";
            }
        }
        qDebug() << "";

        cgsize_t np = numberOfElements(cgio_num, node_id);
        qDebug() << "numberOfElements = " << np;

        std::vector<char> values(data_size+1);
        cgio_read_all_data_type (cgio_num, node_id, data_type, values.data());

        std::string valuestr;

        printNestedFormatValues( data_type, values, ndim, dims, np, valuestr );

        int nchild =-1;
        cgio_number_children (cgio_num, childID, &nchild);

        QStandardItem *childItem = new QStandardItem(childName);
        pItem->appendRow(childItem);

        CgnsNode * node = new CgnsNode();
        node->parent_name = parent_name;
        node->name = childName;
        node->label = label;
        node->item = childItem;
        node->data_type = data_type;
        node->dimstr = dimstr;
        node->data_size = data_size;
        node->valuestr = valuestr;
        node->link_file = link_file;
        node->link_node = link_node;

        nodeMap.insert(std::make_pair(childID, node));
        std::map<double, CgnsNode *>::iterator iter;
        iter = nodeMap.find( childID );
        qDebug() << "iter->second->name = " << iter->second->name;
        qDebug() << "iter->second->parent_name = " << iter->second->parent_name;
        qDebug() << "iter->second->label = " << iter->second->label;

        cgnsNodeMap.insert(std::make_pair(childItem, node));
        std::map<QStandardItem *, CgnsNode *>::iterator cgnsIter;
        cgnsIter = cgnsNodeMap.find( childItem );
        qDebug() << "cgnsIter->second->name = " << cgnsIter->second->name;
        qDebug() << "cgnsIter->second->parent_name = " << iter->second->parent_name;
        qDebug() << "cgnsIter->second->label = " << cgnsIter->second->label;
        qDebug() << "cgnsIter->second->data_type = " << cgnsIter->second->data_type;

        if( nchild <= 0 )
        {
            childItem->setIcon(fileIcon);
        }
        else
        {
            childItem->setIcon(folderIcon);
            if( link_len > 0 )
            {
                childItem->setIcon(forwardIcon);
            }
            ReadCgnsChild( cgio_num, childID, childItem );
        }
    }
}

void CgnsView::ReadCgnsToQStandardItem(const QString &fileName)
{
    int cgio_num = -1;
    cgio_open_file (fileName.toLatin1().data(), CG_MODE_READ, CGIO_FILE_HDF5, &cgio_num);
    ::cgioNum = cgio_num;
    qDebug() << "cgio_num = " << cgio_num;
    double root_id=-1;
    cgio_get_root_id(cgio_num,&root_id);
    qDebug() << "root_id = " << root_id;
    qDebug("root_id = %30.22e", root_id);

    char root_name[33];
    cgio_get_name(cgio_num, root_id, root_name);
    qDebug() << "root_name = " << root_name;

    char root_label[33];
    cgio_get_label (cgio_num, root_id, root_label);
    qDebug() << "root_label = " << root_label;

    int link_len = -1;
    cgio_is_link (cgio_num, root_id, &link_len);
    qDebug() << "link_len = " << link_len;

    std::string link_file;
    std::string link_node;

    if ( link_len > 0 )
    {
        int file_len = -1;
        int name_len = -1;
        cgio_link_size (cgio_num, root_id, &file_len, &name_len);
        qDebug() << "file_len = " << file_len << " name_len = " << name_len;
        std::vector<char> filename_vchar(file_len);
        std::vector<char> name_in_file_vchar(name_len);
        cgio_get_link (cgio_num, root_id, filename_vchar.data(), name_in_file_vchar.data());
        link_file = filename_vchar.data();
        link_node = name_in_file_vchar.data();

        qDebug() << "link_file = " << link_file << " link_node = " << link_node;
    }


    char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
    cgio_get_data_type (cgio_num, root_id, data_type);
    qDebug() << "data_type = " << data_type;

    cglong_t data_size;
    cgio_get_data_size(cgio_num, root_id, &data_size);
    qDebug() << "data_size = " << data_size;

    int ndim;
    cgsize_t dims[CGIO_MAX_DIMENSIONS];
    cgio_get_dimensions (cgio_num, root_id, &ndim, dims);
    qDebug() << "ndim = " << ndim;
    std::string dimstr;

    qDebug() << "dims = ";
    for( int i = 0; i < ndim; ++ i )
    {
        qDebug() << dims[i];
        dimstr += std::to_string(dims[i]);
        if ( i != ndim-1 )
        {
            dimstr += " ";
        }
    }
    qDebug() << "";

    cgsize_t np = numberOfElements(cgio_num, root_id);
    qDebug() << "numberOfElements = " << np;

    std::vector<char> values(data_size+1);
    cgio_read_all_data_type (cgio_num, root_id, data_type, values.data());
    std::string valuestr;

    printNestedFormatValues( data_type, values, ndim, dims, np, valuestr );

    char version[33];
    sprintf (version, "%g", CGNS_DOTVERS);
    qDebug() << "version = " << version;

    char lib_version[CGIO_MAX_VERSION_LENGTH+1];
    cgio_library_version (cgio_num, lib_version);
    qDebug() << "lib_version = " << lib_version;

    QIcon folderIcon(":/folder.png");
    QIcon fileIcon(":/file.jpg");

    QStandardItemModel *itemModel = static_cast<QStandardItemModel *>(this->treeView->model());
    QStandardItem *rootItem = new QStandardItem( "/" );
    itemModel->appendRow(rootItem);
    rootItem->setIcon(folderIcon);

    CgnsNode * node = new CgnsNode();
    node->parent_name = "";
    node->name = "/";
    node->label = root_label;
    node->item = rootItem;
    node->data_type = data_type;
    node->dimstr = dimstr;
    node->data_size = data_size;
    node->valuestr = valuestr;
    node->link_file = link_file;
    node->link_node = link_node;

    nodeMap.insert(std::make_pair(root_id, node));
    std::map<double, CgnsNode *>::iterator iter;
    iter = nodeMap.find( root_id );
    qDebug() << "iter->second->name = " << iter->second->name;
    qDebug() << "iter->second->parent_name = " << iter->second->parent_name;
    qDebug() << "iter->second->label = " << iter->second->label;

    cgnsNodeMap.insert(std::make_pair(rootItem, node));
    std::map<QStandardItem *, CgnsNode *>::iterator cgnsIter;
    cgnsIter = cgnsNodeMap.find( rootItem );
    qDebug() << "cgnsIter->second->name = " << cgnsIter->second->name;
    qDebug() << "cgnsIter->second->parent_name = " << cgnsIter->second->parent_name;
    qDebug() << "cgnsIter->second->label = " << cgnsIter->second->label;
    qDebug() << "cgnsIter->second->data_type = " << cgnsIter->second->data_type;

    ReadCgnsChild(cgio_num,root_id,rootItem);
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
    rootItem->appendRow(libraryVersionItem);
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
        rootItem->appendRow(baseItem);
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
        qDebug() <<  "bcItem->index()=" << bcItem->index();
        qDebug() <<  "bcItem=" << bcItem;
        qDebug() <<  "bcItem->model()->itemFromIndex(bcItem->index())=" << bcItem->model()->itemFromIndex(bcItem->index());
        qDebug() <<  "bcItem->index().internalPointer()=" << bcItem->index().internalPointer();
        qDebug() <<  "bcItem->model()=" << bcItem->model();
        qDebug() <<  "bcItem->parent()=" << bcItem->parent();
        qDebug() <<  "zoneBcItem=" << zoneBcItem;

        QStandardItem *pointRangeItem = new QStandardItem( "PointRange" );
        bcItem->appendRow(pointRangeItem);
        pointRangeItem->setIcon(fileIcon);

        qDebug() <<  "pointRangeItem->index()=" << pointRangeItem->index();
        qDebug() <<  "pointRangeItem=" << pointRangeItem;

        QStandardItem *gridLocationItem = new QStandardItem( "GridLocation" );
        bcItem->appendRow(gridLocationItem);
        gridLocationItem->setIcon(fileIcon);
    }
}

void CgnsView::SetCgnsPanel(CgnsPanel * cgnsPanel)
{
    this->cgnsPanel = cgnsPanel;
}

void CgnsView::onTreeItemClicked(const QModelIndex &index)
{
    qDebug() <<"CgnsView::on_treeView_clicked" << "index=" << index;
    int row = index.row();
    int col = index.column();
    qDebug() <<"item row=" << row << "column=" << col;

    QString text = index.data(Qt::DisplayRole).toString();
    //this->cgnsPanel->Display(text);
    qDebug() << "text = " << text;

    QStandardItemModel * itemModel = static_cast<QStandardItemModel *>( this->treeView->model() );
    QStandardItem * item = itemModel->itemFromIndex( index );

    std::map<QStandardItem *, CgnsNode *>::iterator cgnsIter;
    cgnsIter = cgnsNodeMap.find( item );
    qDebug() << "cgnsIter->second->name = " << cgnsIter->second->name;
    qDebug() << "cgnsIter->second->parent_name = " << cgnsIter->second->parent_name;
    qDebug() << "cgnsIter->second->label = " << cgnsIter->second->label;

    this->cgnsPanel->DisplayNode( cgnsIter->second );


    // this->cgnsPanel->Display(node);
    //DisplayNode(CgnsNode *node)

}

