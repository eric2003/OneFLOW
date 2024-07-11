#include "Terminal.h"
#include "./ui_terminal.h"
#include <QKeyEvent>
#include <QDir>
#include <QTextDocument>
#include <QTextBlock>
#include <QScrollBar>
#include <QCompleter>

Terminal::Terminal(QWidget *parent)
    : QWidget{parent}
    , ui(new Ui::Terminal)
{
    ui->setupUi(this);
    //currentPath = QCoreApplication::applicationDirPath();
    currentPath =  QDir::currentPath();
    qDebug() << "currentPath="<<currentPath;

    wordList << "apple" << "banana" << "cherry" << "date";

    this->ui->textEdit->append( "1" );
    this->ui->textEdit->setPlainText( wordList[0] );
    this->ui->textEdit->append( "2" );
    this->ui->textEdit->setPlainText( wordList[1] );


    QTextCursor cursor = this->ui->textEdit->textCursor();
    int cursorPosition = cursor.position();
    this->ui->textEdit->setCursorWidth(4);
    this->ui->textEdit->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    this->ui->textEdit->moveCursor(QTextCursor::End);
    cursor.movePosition(QTextCursor::Right);

    int lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    qDebug() << "cursorPosition="<<cursorPosition;
    qDebug() << "lineNumber="<<lineNumber;

    // foreach(QString str, wordList)
    // {
    //     qDebug() << str;
    //     this->ui->textEdit->append( str );
    // }

    QString currStr = QDir::currentPath() + ">";
    this->lastCommand = currStr;
    QString currStr1 = QDir::toNativeSeparators(currStr);
    this->ui->textEdit->append( currStr );
    cursorPosition = cursor.position();
    lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    qDebug() << "cursorPosition="<<cursorPosition;
    qDebug() << "lineNumber="<<lineNumber;

    this->ui->textEdit->append( currStr1 );
    this->lastCommand = currStr1;

    cursorPosition = cursor.position();
    lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    qDebug() << "cursorPosition="<<cursorPosition;
    qDebug() << "lineNumber="<<lineNumber;

    QPalette p = this->ui->textEdit->palette(); // define pallete for textEdit..
    QColor bkcolor(12, 12, 12);
    //QColor bkcolor(255, 255, 255);
    QColor tcolor(204, 204, 204);
    //p.setColor(QPalette::Base, Qt::black); // set color "Red" for textedit base
    //p.setColor(QPalette::Text, Qt::white); // set text color which is selected from color pallete
    p.setColor(QPalette::Base, bkcolor);
    p.setColor(QPalette::Text, tcolor);
    this->ui->textEdit->setPalette(p);
    QFont font("NSimSun", 12);
    this->ui->textEdit->setFont(font);

    ui->textEdit->installEventFilter(this);
    ui->textEdit->viewport()->installEventFilter(this);

    QScrollBar * vBar = this->ui->textEdit->verticalScrollBar();
    //vBar->setStyleSheet("QScrollBar:vertical { background-color:rgb(240,240,240);width: 10px; }");
    vBar->setStyleSheet("QSlider::groove { background: transparent; width: 10px; } ");
    this->procCmd = new QProcess( this );
    //Command line related
    QObject::connect( this->procCmd, &QProcess::readyRead, this, &Terminal::ReadOutput );                //Read command line data
    QObject::connect( this->procCmd, &QProcess::readyReadStandardOutput, this, &Terminal::ReadStandardOutput );  //Read command line standard data
    QObject::connect( this->procCmd, &QProcess::finished, this,  &Terminal::FinishedProcess );           //Command line processing ends
    QObject::connect( this->procCmd, &QProcess::errorOccurred, this, &Terminal::ErrorProcess );          //Command line error handling
    QObject::connect( this->procCmd, &QProcess::stateChanged, this, &Terminal::OnStateChanged );          //Command line error handling

    //this->procCmd->start("cmd", QStringList() << "/c" << "ping www.baidu.com" );
}

bool Terminal::eventFilter(QObject *obj, QEvent *event)
{
    //qDebug() << "Terminal event->type()=" << event->type();
    if (event->type() == QEvent::KeyPress)
    {
        QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
        if (keyEvent->key() == Qt::Key_Return)
        {
            qDebug() << "Terminal Enter key pressed!";
            // 执行回车键按下后的操作
            this->OnReturnKeyPressed();
            return true;
        }
        else if( keyEvent->key() == Qt::Key_Backspace ||
                 keyEvent->key() == Qt::Key_Left )
        {
            qDebug() << "eventFilter lastCommand=" << lastCommand;
            this->currentPath =  QDir::currentPath();
            qDebug() << "eventFilter currentPath=" << QDir::currentPath();
            QTextCursor cursor = this->ui->textEdit->textCursor();
            this->OnKeyBackPressed();
            if( lastCommand.length() <= currentPath.length()+1)
            {
                return true;
            }
        }
        else if( keyEvent->key() == Qt::Key_Up )
        {
            int N = wordList.size();
            QString s = wordList[(word_index++)%N];

            QString currStr = QDir::currentPath() + ">";

            QTextCursor cursor = this->ui->textEdit->textCursor();
            int cursorPosition = cursor.position();
            qDebug() << "cursorPosition=" << cursor.position();
            cursor.movePosition( QTextCursor::StartOfLine );
            int startPosition = cursor.position() + currStr.length();
            qDebug() << "cursorPosition=" << cursor.position();
            cursor.setPosition( startPosition );
            cursor.movePosition( QTextCursor::EndOfLine, QTextCursor::KeepAnchor );
            cursor.deleteChar();
            cursor.insertText( s );
            this->ui->textEdit->setTextCursor(cursor);
            return true;
        }
        else if (keyEvent->key() == Qt::Key_Tab)
        {
            qDebug() << "Terminal Tab key pressed!";
            this->ProcessTabKey();
            return true;
        }
    }
    else if( event->type() == QEvent::MouseButtonPress )
    {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
        qDebug() << "eventFilter QEvent::MouseButtonPress";
        qDebug() << "eventFilter mouseEvent->pos()="<<mouseEvent->pos();
        QTextCursor mouseCursor = this->ui->textEdit->cursorForPosition(mouseEvent->pos());
        int mouseCursorPosition = mouseCursor.position();
        int mouseCursorLineNumber = this->ui->textEdit->document()->findBlock(mouseCursorPosition).blockNumber() + 1;
        qDebug() << "QEvent::MouseButtonPress mouseCursorPosition="<<mouseCursorPosition;
        qDebug() << "QEvent::MouseButtonPress mouseCursorLineNumber="<<mouseCursorLineNumber;

        QTextCursor currentCursor = this->ui->textEdit->textCursor();
        int currentCursorPosition = currentCursor.position();
        int currentCursorLineNumber = this->ui->textEdit->document()->findBlock(currentCursorPosition).blockNumber() + 1;

        qDebug() << "QEvent::MouseButtonPress currentCursorPosition="<<currentCursorPosition;
        qDebug() << "QEvent::MouseButtonPress currentCursorLineNumber="<<currentCursorLineNumber;

        this->ui->textEdit->moveCursor(QTextCursor::End);
        QTextCursor currentCursor1 = this->ui->textEdit->textCursor();
        int currentCursorPosition1 = currentCursor1.position();
        int currentCursorLineNumber1 = this->ui->textEdit->document()->findBlock(currentCursorPosition1).blockNumber() + 1;

        qDebug() << "QEvent::MouseButtonPress currentCursorPosition1="<<currentCursorPosition1;
        qDebug() << "QEvent::MouseButtonPress currentCursorLineNumber1="<<currentCursorLineNumber1;

        if ( mouseEvent->button() == Qt::LeftButton )
        {
            qDebug() << "eventFilter Qt::LeftButton";
            QWidget::mousePressEvent(mouseEvent);
        }
        else if ( mouseEvent->button() == Qt::RightButton )
        {
            qDebug() << "eventFilter Qt::RightButton";
        }
    }
    else if( event->type() == QEvent::MouseButtonRelease )
    {
        qDebug() << "eventFilter QEvent::MouseButtonRelease";
    }

    return QObject::eventFilter(obj, event);
}

void Terminal::resizeEvent(QResizeEvent *event)
{
    qDebug() << "Terminal::resizeEvent(QResizeEvent *event)";
    qDebug() << "Terminal::resizeEvent this->rect()=" << this->rect();
    qDebug() << "Terminal::resizeEvent event->size()=" << event->size();
    int borderWidth = 10;

    this->ui->textEdit->resize(event->size().width()-borderWidth, event->size().height()-borderWidth);
}

void Terminal::OnReturnKeyPressed()
{
    QTextCursor cursor = this->ui->textEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    int cursorPosition = cursor.position();
    QTextDocument *textDocument = this->ui->textEdit->document();

    QTextBlock textBlock = textDocument->findBlock(cursorPosition);
    QString selectLine = textBlock.text();
    this->lastCommand = selectLine;
    qDebug() << "Terminal::OnReturnKeyPressed selectLine="<<selectLine;
    qDebug() << "Terminal::OnReturnKeyPressed this->lastCommand="<<this->lastCommand;

    QString currStr = QDir::currentPath() + ">";

    QStringList argument;
    QString cmdStr = selectLine.mid(currStr.length());

    qDebug() << "Terminal::OnReturnKeyPressed cmdStr="<<cmdStr;
    qDebug() << "Terminal::OnReturnKeyPressed selectLine="<<selectLine;

    //QString cmd = "cd d:\\ && dir";
    Analysis( cmdStr );
    argument << "/c" << cmdStr;
    this->procCmd->start("cmd", argument );
    if ( !this->procCmd->waitForStarted() ) {
        qDebug() << "Terminal::OnReturnKeyPressed waitForStarted Failed";
        return;
    }

    if ( !this->procCmd->waitForFinished() ) {
        qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished Failed";
        return;
    }

    currStr = QDir::currentPath() + ">";
    this->ui->textEdit->append( currStr );
    this->lastCommand = currStr + cmdStr;

    //Information Output
    qDebug() << "Success:OnReturnKeyPressed";
}

void Terminal::OnKeyBackPressed()
{
    QTextCursor cursor = this->ui->textEdit->textCursor();
    QTextBlock textBlock = cursor.block();

    QString selectLine = textBlock.text();
    this->lastCommand = selectLine;
}

QString Terminal::GetCommandString()
{
    QTextCursor cursor = this->ui->textEdit->textCursor();
    QTextBlock textBlock = cursor.block();

    QString selectLine = textBlock.text();
    QString currStr = QDir::currentPath() + ">";
    QString cmdStr = selectLine.mid(currStr.length());
    QString trimedCmd = cmdStr.trimmed();
    qDebug() << "cmdStr="<<cmdStr;
    qDebug() << "trimedCmd="<<trimedCmd;
    return trimedCmd;
}

void Terminal::CommandCompletion( const QString &str )
{
    QTextCursor cursor = this->ui->textEdit->textCursor();
    cursor.insertText( str );
    this->ui->textEdit->setTextCursor(cursor);
}

QString Terminal::FindCommonString(const QString &stra, const QString &strb)
{
    int lena = stra.length();
    int lenb = strb.length();
    int len = std::min(lena, lenb);
    qDebug() << "len=" << len;
    qDebug() << "lena=" << lena;
    qDebug() << "lenb=" << lenb;
    qDebug() << "stra=" << stra;
    qDebug() << "strb=" << strb;
    int ipos = -1;
    for ( int i = 0; i < len; ++ i )
    {
        if ( stra[i] != strb[i] )
        {
            break;
        }
        ipos = i;
    }
    qDebug() << "ipos=" << ipos;
    QString cs = "";
    if( ipos != -1 )
    {
        cs = stra.mid(0, ipos+1);
    }
    return cs;
}

bool Terminal::InStringList(QChar a, int ipos, const QStringList &strlist)
{
    for( int i = 0; i < strlist.size(); ++ i )
    {
        if( a != strlist[i][ipos] ) return false;
    }
    return true;
}

QString Terminal::FindCommonString(const QStringList &strlist)
{
    int len = std::numeric_limits<int>::max();
    for( int i = 0; i < strlist.size(); ++ i )
    {
        len = std::min(len, static_cast<int>(strlist[i].length()));
    }
    qDebug() << "len=" << len;

    const QString & stra = strlist[0];
    int ipos = -1;
    for( int i = 0; i < len; ++ i )
    {
        bool flag = InStringList(stra[i], i, strlist);
        if ( ! flag ) break;
        ipos = i;
    }
    QString cs = "";
    if ( ipos != -1 )
    {
        cs = stra.mid(0,ipos+1);
    }
    return cs;
}

void Terminal::ProcessTabKey()
{
    QString cmdStr = this->GetCommandString();
    qDebug() << "Terminal::ProcessTabKey cmdStr="<<cmdStr;

    QDir dir;
    dir.setFilter( QDir::Files | QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot );
    QFileInfoList filelist = dir.entryInfoList();
    qDebug() << filelist;

    QString str = "Hello world!";
    int index = str.indexOf("world");
    qDebug() << "Hello world! str.indexOf(\"world\") index=" << index;
    int index1 = str.indexOf("world");
    QString c1 = str.mid(0,2);
    qDebug() << "c1=" << c1;
    QStringList cmdlist = cmdStr.split(" ");
    qDebug() << "cmdlist = " << cmdlist;
    QString cmdPara = "";
    if( cmdlist.size() >= 2 )
    {
        cmdPara = cmdlist.at(1);
    }

    qDebug() << "cmdPara = " << cmdPara;
    QStringList candidateList;
    QStringList candidateResList;

    for ( QFileInfo & fileInfo : filelist )
    {
        qDebug() << "filename=" << fileInfo.fileName();
        QString fileName = fileInfo.fileName();
        int len = cmdPara.length();
        QString res = fileName.mid(0,len);
        QString res1 = fileName.mid(len);
        qDebug() << "res" << res << " res1 = " << res1;
        if ( cmdPara == res )
        {
            qDebug() << "cmdPara == res fileName = " << fileName << " res = " << res;
            candidateList.append( fileName );
            candidateResList.append( res1 );
        }
        else
        {
            qDebug() << "cmdPara /= res fileName = " << fileName << " res = " << res;
        }
    }
    qDebug() << "candidateList = " << candidateList;
    QString co = FindCommonString("test", "ttt");
    QString co1 = FindCommonString(candidateList);
    QString co2 = FindCommonString(candidateResList);
    qDebug() << " co = " << co;
    qDebug() << " co1 = " << co1;
    qDebug() << " co2 = " << co2;

    if ( !co2.isEmpty())
    {
        this->CommandCompletion(co2);
    }

    if ( candidateList.size() == 1 )
    {
        //this->CommandCompletion(candidateResList[0]);
    }
}

void Terminal::Analysis( QString & cmdString )
{
    qDebug() << "cmdString="<< cmdString;
    QString c1 = cmdString.mid(0,3);
    qDebug() << "c1="<< c1;
    if ( cmdString.mid(0,3) == "cd " )
    {
        QString targetDir = cmdString.mid(3);
        qDebug() << "cd command";
        qDebug() << "targetDir="<< targetDir;
        QString currentPath1 = QDir::currentPath();
        qDebug() << "currentPath1="<< currentPath1;
        QDir::setCurrent(targetDir);
        QString currentPath2 = QDir::currentPath();
        qDebug() << "currentPath2="<< currentPath2;

    }
}

void Terminal::OnStateChanged()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Success:OnStateChanged";
}

void Terminal::ReadOutput()
{
    qDebug() << Q_FUNC_INFO;

    //QByteArray qByteRead = this->procCmd->readAll() + this->procCmd->readAllStandardOutput();
    QByteArray qByteRead = this->procCmd->readAllStandardOutput();

    if ( qByteRead.isEmpty() ) return;

    QString word = QString::fromLocal8Bit( qByteRead );
    int pos = word.lastIndexOf("\r\n");
    if (pos != -1) {
        word.remove(pos, 2); // 移除末尾的 "\r\n"
    }
    //qDebug() << "word="<< word;
    this->ui->textEdit->append( word );

    QTextCursor cursor = this->ui->textEdit->textCursor();
    int cursorPosition = cursor.position();
    int lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    //qDebug() << "cursorPosition="<<cursorPosition;
    //qDebug() << "lineNumber="<<lineNumber;
}

void Terminal::ReadStandardOutput()
{
    qDebug() << Q_FUNC_INFO;

    QByteArray qByteRead = this->procCmd->readAllStandardOutput();
    qDebug() << qByteRead;
    if ( qByteRead.isEmpty() ) return;

    this->ui->textEdit->append( QString::fromLocal8Bit( qByteRead ).trimmed() );
}

void Terminal::FinishedProcess()
{
    //Receive data
    int flag = this->procCmd->exitCode();

    //Information Output
    qDebug() << "Success:FinishedProcess():" << flag;
}

void Terminal::ErrorProcess()
{
    //Receive data
    int err_code  = this->procCmd->exitCode();
    QString err = this->procCmd->errorString();

    //Display Data
    this->ui->textEdit->append(QString("error code:%1").arg(err_code));
    this->ui->textEdit->append(err);

    //Information Output
    qDebug() << "Success:ErrorProcess():" << err;
}

