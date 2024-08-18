#include "explorer.h"
#include <QTreeView>
#include <QStandardItemModel>
#include <QHeaderView>
#include <QDir>
#include <QDirIterator>
#include <QResizeEvent>

Explorer::Explorer(QWidget *parent)
    : QWidget{parent}
{
    this->treeView = new QTreeView(this);
    this->treeView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    //this->treeView->setGeometry(parent->rect());

    // this->treeView->setStyleSheet("QTreeView::branch:has-children:!has-siblings:closed,\
    // QTreeView::branch:closed:has-children:has-siblings{border-image: none; image: url(:/plus.png);}\
    // QTreeView::branch:open:has-children:!has-siblings,\
    // QTreeView::branch:open:has-children:has-siblings{border-image: none; image: url(:/minus2.png);}");

    this->SetUpExplorerInfo();
}

void Explorer::resizeEvent(QResizeEvent *event)
{
    this->treeView->resize(event->size());

}

void Explorer::SetUpHeader()
{
    QStandardItemModel *itemModel = new QStandardItemModel();
    this->treeView->setModel(itemModel);
    this->treeView->header()->setSectionResizeMode(QHeaderView::Stretch);
    this->treeView->show();
    QStringList labels;
    labels << "EXPLORER";
    itemModel->setHorizontalHeaderLabels(labels);
}

void Explorer::SetUpExplorerInfo()
{
    this->SetUpHeader();
    //PrintDirInfo();

    QStandardItemModel *itemModel = static_cast<QStandardItemModel *>(this->treeView->model());

    QStandardItem *rootItem = new QStandardItem(QDir::currentPath());
    itemModel->appendRow(rootItem);

    this->treeView->expandAll();

    AddExplorerItem(rootItem, ".");
}

void Explorer::PrintDirInfo()
{
    QDirIterator dirIterator(".", QDir::Files | QDir::NoDotAndDotDot, QDirIterator::Subdirectories);
    while (dirIterator.hasNext())
    {
        dirIterator.next(); // 获取当前路径
        QFileInfo fileInfo = dirIterator.fileInfo();
        qDebug() << "文件路径：" << fileInfo.filePath();
        qDebug() << "文件大小：" << fileInfo.size();
        qDebug() << "是否是目录：" << fileInfo.isDir();
        qDebug() << "fileName()" << fileInfo.fileName();
    }
}

void Explorer::AddExplorerItem(QStandardItem *rootItem)
{
    QDir dir = QDir::current();
    qDebug() << "dir="<<dir;

    qDebug() << "dir.dirName()="<<dir.dirName();

    dir.setFilter(QDir::Dirs|QDir::Files);

    QFileInfoList fileInfoLists = dir.entryInfoList();

    for( int i = 0; i < fileInfoLists.size(); ++ i )
    {
        qDebug() << "fileName:" << i << " = " << fileInfoLists[i].fileName();
        qDebug() << "baseName:" << i << " = " << fileInfoLists[i].baseName();
        qDebug() << "completeBaseName:" << i << " = " << fileInfoLists[i].completeBaseName();
        qDebug() << "filePath:" << i << " = " << fileInfoLists[i].filePath();
        qDebug() << "absoluteFilePath:" << i << " = " << fileInfoLists[i].absoluteFilePath();
        qDebug() << "canonicalFilePath:" << i << " = " << fileInfoLists[i].canonicalFilePath();
        qDebug() << "filesystemFilePath:" << i << " = " << fileInfoLists[i].filesystemFilePath();


        bool fileflag = fileInfoLists[i].isFile();
        bool dirflag = fileInfoLists[i].isDir();
        QString fn = fileInfoLists[i].fileName();

        if ( dirflag )
        {
            fn = "+" + fn;
        }
        QStandardItem *childItem = new QStandardItem(fn);
        rootItem->appendRow(childItem);
    }
}

void Explorer::AddExplorerItem(QStandardItem *rootItem, const QString &filename)
{
    QDir dir(filename);
    //qDebug() << "AddExplorerItem(QStandardItem *rootItem, const QString &filename) dir="<<dir << "filename="<<filename;
    //qDebug() << "QDir::currentPath()=" << QDir::currentPath();

    //qDebug() << "dir.dirName()="<<dir.dirName();

    dir.setFilter(QDir::Dirs|QDir::Files| QDir::NoDotAndDotDot);

    QFileInfoList fileInfoLists = dir.entryInfoList();
    //qDebug() << "fileInfoLists.size(): " << fileInfoLists.size();

    for( int i = 0; i < fileInfoLists.size(); ++ i )
    {
        //qDebug() << "fileName:" << i << " = " << fileInfoLists[i].fileName();
        //qDebug() << "baseName:" << i << " = " << fileInfoLists[i].baseName();
        //qDebug() << "suffix:" << i << " = " << fileInfoLists[i].suffix();
        //qDebug() << "completeBaseName:" << i << " = " << fileInfoLists[i].completeBaseName();
        //qDebug() << "filePath:" << i << " = " << fileInfoLists[i].filePath();
        //qDebug() << "absoluteFilePath:" << i << " = " << fileInfoLists[i].absoluteFilePath();
        //qDebug() << "canonicalFilePath:" << i << " = " << fileInfoLists[i].canonicalFilePath();
        //qDebug() << "filesystemFilePath:" << i << " = " << fileInfoLists[i].filesystemFilePath();


        bool fileflag = fileInfoLists[i].isFile();
        bool dirflag = fileInfoLists[i].isDir();
        //qDebug() << "dirflag: " << dirflag;
        //qDebug() << "fileflag: " << fileflag;
        QString fn = fileInfoLists[i].fileName();

        QStandardItem *childItem = new QStandardItem(fn);

        if ( dirflag )
        {
            childItem->setIcon(QIcon(":/folder.png"));
        }
        else
        {
            QString suf = fileInfoLists[i].suffix();
            if( suf ==  "cmake")
            {
                childItem->setIcon(QIcon(":/cmake.jpg"));
            }
            else
            {
                childItem->setIcon(QIcon(":/file.jpg"));
            }
        }

        rootItem->appendRow(childItem);

        if ( dirflag )
        {
            this->AddExplorerItem(childItem, fileInfoLists[i].absoluteFilePath());
        }
    }
}

