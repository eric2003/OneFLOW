#include "Terminal.h"
#include "./ui_terminal.h"
#include <QKeyEvent>
#include <QDir>
#include <QTextEdit>
#include <QTextDocument>
#include <QTextBlock>
#include <QScrollBar>
#include <QCompleter>

QTextEdit *textEdit;

Terminal::Terminal(QWidget *parent)
    : QWidget{parent}
    , ui(new Ui::Terminal)
{
    ui->setupUi(this);

    ::textEdit = ui->textEdit;

    QTextCursor cursor = this->ui->textEdit->textCursor();
    int cursorPosition = cursor.position();
    this->ui->textEdit->setCursorWidth(4);
    this->ui->textEdit->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    this->ui->textEdit->moveCursor(QTextCursor::End);
    cursor.movePosition(QTextCursor::Right);

    int lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    qDebug() << "cursorPosition="<<cursorPosition;
    qDebug() << "lineNumber="<<lineNumber;

    QString currStr = QDir::currentPath() + ">";
    this->ui->textEdit->append( currStr );
    //this->terminalThread->lastCommand = currStr;

    cursorPosition = cursor.position();
    lineNumber = this->ui->textEdit->document()->findBlock(cursorPosition).blockNumber() + 1;
    qDebug() << "cursorPosition="<<cursorPosition;
    qDebug() << "lineNumber="<<lineNumber;

    QPalette p = this->ui->textEdit->palette(); // define pallete for textEdit..
    //QColor bkcolor(12, 12, 12);
    //QColor tcolor(204, 204, 204);
    QColor bkcolor(255, 255, 255);
    QColor tcolor(12, 12, 12);
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
    vBar->setStyleSheet("QSlider::groove { background: transparent; width: 10px; } ");

    //this->procCmd->start("cmd", QStringList() << "/c" << "ping www.baidu.com" );

    this->procCmd = new QProcess();

    //Command line related
    QObject::connect( this->procCmd, &QProcess::started, this, &Terminal::OnStarted );
    QObject::connect( this->procCmd, &QProcess::readyRead, this, &Terminal::ReadOutput );                //Read command line data
    QObject::connect( this->procCmd, &QProcess::readyReadStandardOutput, this, &Terminal::ReadStandardOutput );  //Read command line standard data
    QObject::connect( this->procCmd, &QProcess::finished, this,  &Terminal::FinishedProcess );           //Command line processing ends
    QObject::connect( this->procCmd, &QProcess::errorOccurred, this, &Terminal::ErrorProcess );          //Command line error handling
    QObject::connect( this->procCmd, &QProcess::stateChanged, this, &Terminal::OnStateChanged );
}

bool Terminal::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::KeyPress)
    {
        QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
        // this->terminalThread->flagOnReturnKeyPressed = false;
        // this->terminalThread->flagOnKeyBackPressed = false;
        // this->terminalThread->flagOnKeyUpPressed = false;
        // this->terminalThread->flagProcessTabKey = false;
        if (keyEvent->key() == Qt::Key_Return)
        {
            qDebug() << "Terminal Enter key pressed!";
            //this->terminalThread->OnReturnKeyPressed();
            this->OnReturnKeyPressed();
            //this->terminalThread->flagOnReturnKeyPressed = true;
            //this->terminalThread->start();
            return true;
        }
        else if( keyEvent->key() == Qt::Key_Backspace ||
            keyEvent->key() == Qt::Key_Left )
        {
            qDebug() << "eventFilter currentPath=" << QDir::currentPath();
            QTextCursor cursor = this->ui->textEdit->textCursor();
            this->OnKeyBackPressed();
            if( this->lastCommand.length() <= QDir::currentPath().length()+1)
            {
                return true;
            }
        }
        else if( keyEvent->key() == Qt::Key_Up )
        {
            this->OnKeyUpPressed();
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

        if ( mouseEvent->button() == Qt::LeftButton )
        {
            qDebug() << "eventFilter Qt::LeftButton";
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

    qDebug() << "Terminal::resizeEvent this->width()=" << this->width();
    qDebug() << "Terminal::resizeEvent this->height()=" << this->height();

    int borderWidth = 0;
    int yTop = 30;

    int layoutWidth = this->width() - borderWidth;
    int layoutHeight = this->height() - yTop - borderWidth;

    this->ui->verticalLayoutWidget->setGeometry( QRect(0, yTop, layoutWidth, layoutHeight) );
}

void Terminal::OnStarted()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Terminal::OnStarted()+++++++++++++++++++++++++++++++++++";
}

void Terminal::ReadOutput()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Terminal::ReadOutput()+++++++++++++++++++++++++++++++++++";
    QByteArray qByteRead = this->procCmd->readAllStandardOutput();

    if ( qByteRead.isEmpty() ) return;

    QString word = QString::fromLocal8Bit( qByteRead );
    int pos = word.lastIndexOf("\r\n");
    if (pos != -1) {
        word.remove(pos, 2); // 移除末尾的 "\r\n"
    }
    qDebug() << "ReadOutput word="<< word;

    ::textEdit->append( word );
}

void Terminal::ReadStandardOutput()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "MainWindow::ReadStandardOutput()+++++++++++++++++++++++++++++++++++";
    QByteArray qByteRead = this->procCmd->readAllStandardOutput();
    qDebug() << qByteRead;
    if ( qByteRead.isEmpty() ) return;
    qDebug() << "ReadStandardOutput qByteRead="<< qByteRead;

    ::textEdit->append( QString::fromLocal8Bit( qByteRead ).trimmed() );
}

void Terminal::FinishedProcess()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Terminal::FinishedProcess()+++++++++++++++++++++++++++++++++++";

    //Receive data
    int flag = this->procCmd->exitCode();

    //Information Output
    qDebug() << "Success:FinishedProcess(): this->procCmd->exitCode() = " << flag;

    // qDebug() << "Terminal::OnReturnKeyPressed waitForStarted";
    // if ( !this->procCmd->waitForStarted() ) {
    //     qDebug() << "Terminal::OnReturnKeyPressed waitForStarted Failed";
    //     return;
    // }

    // qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished Begin";

    // if ( !this->procCmd->waitForFinished() ) {
    //     qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished Failed";
    //     return;
    // }

    //qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished End";

    QString  currStr = QDir::currentPath() + ">";
    qDebug() << "currStr=" << currStr;
    qDebug() << "this->cmdStr=" << this->cmdStr;
    qDebug() << "::textEdit=" << ::textEdit;

    ::textEdit->append( currStr );
    this->lastCommand = currStr + this->cmdStr;
    qDebug() << "this->lastCommand=" << this->lastCommand;

    if ( !cmdStr.isEmpty() )
    {
        qDebug() << "enter !cmdStr.isEmpty()";
        this->historyCmdList.append( this->cmdStr );
        qDebug() << "this->historyCmdList = " << this->historyCmdList;
    }
}

void Terminal::ErrorProcess()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Terminal::ErrorProcess()+++++++++++++++++++++++++++++++++++";

    //Receive data
    int err_code  = this->procCmd->exitCode();
    QString err = this->procCmd->errorString();

    //Display Data
    ::textEdit->append(QString("error code:%1").arg(err_code));
    ::textEdit->append(err);

    //Information Output
    qDebug() << "Success:ErrorProcess():" << err;
}

void Terminal::OnStateChanged()
{
    qDebug() << Q_FUNC_INFO;
    qDebug() << "Terminal::OnStateChanged()+++++++++++++++++++++++++++++++++++";
}

void Terminal::OnReturnKeyPressed()
{
    QTextCursor cursor = ::textEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    int cursorPosition = cursor.position();
    QTextDocument *textDocument = ::textEdit->document();

    QTextBlock textBlock = textDocument->findBlock(cursorPosition);
    QString selectLine = textBlock.text();
    this->lastCommand = selectLine;
    qDebug() << "Terminal::OnReturnKeyPressed selectLine="<<selectLine;
    qDebug() << "Terminal::OnReturnKeyPressed this->lastCommand="<<this->lastCommand;

    QString currStr = QDir::currentPath() + ">";

    QStringList argument;
    this->cmdStr = selectLine.mid(currStr.length());

    qDebug() << "Terminal::OnReturnKeyPressed cmdStr="<<this->cmdStr;
    qDebug() << "Terminal::OnReturnKeyPressed selectLine="<<selectLine;

    qDebug() << "Terminal::Analysis Begin";
    Analysis( this->cmdStr );
    qDebug() << "Terminal::Analysis End";
    argument << "/c" << this->cmdStr;
    qDebug() << "argument =" << argument;
    this->procCmd->start("cmd", argument );
    // qDebug() << "Terminal::OnReturnKeyPressed waitForStarted";
    // if ( !this->procCmd->waitForStarted() ) {
    //     qDebug() << "Terminal::OnReturnKeyPressed waitForStarted Failed";
    //     return;
    // }

    // qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished Begin";

    // if ( !this->procCmd->waitForFinished() ) {
    //     qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished Failed";
    //     return;
    // }

    // qDebug() << "Terminal::OnReturnKeyPressed External waitForFinished End";

    // currStr = QDir::currentPath() + ">";
    // qDebug() << "currStr=" << currStr;
    // qDebug() << "cmdStr=" << cmdStr;
    // qDebug() << "this->textEdit=" << this->textEdit;

    // this->textEdit->append( currStr );
    // this->lastCommand = currStr + cmdStr;
    // qDebug() << "this->lastCommand=" << this->lastCommand;

    // if ( !cmdStr.isEmpty() )
    // {
    //     qDebug() << "enter !cmdStr.isEmpty()";
    //     this->historyCmdList.append( cmdStr );
    //     qDebug() << "this->historyCmdList = " << this->historyCmdList;
    // }

    //Information Output
    qDebug() << "Success:OnReturnKeyPressed";
}

void Terminal::OnKeyBackPressed()
{
    QTextCursor cursor = ::textEdit->textCursor();
    QTextBlock textBlock = cursor.block();

    QString selectLine = textBlock.text();
    this->lastCommand = selectLine;
}

void Terminal::OnKeyUpPressed()
{
    int N = this->historyCmdList.size();
    if( N == 0 ) return;

    QString s = this->historyCmdList[(word_index++)%N];

    QString currStr = QDir::currentPath() + ">";

    QTextCursor cursor = ::textEdit->textCursor();
    cursor.movePosition( QTextCursor::StartOfLine );
    int startPosition = cursor.position() + currStr.length();
    qDebug() << "cursorPosition=" << cursor.position();
    cursor.setPosition( startPosition );
    cursor.movePosition( QTextCursor::EndOfLine, QTextCursor::KeepAnchor );
    cursor.deleteChar();
    cursor.insertText( s );
    ::textEdit->setTextCursor(cursor);
}

void Terminal::ProcessTabKey()
{
    QString cmdStr = this->GetCommandString();
    qDebug() << "Terminal::ProcessTabKey cmdStr="<<cmdStr;

    QDir dir;
    dir.setFilter( QDir::Files | QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot );
    QFileInfoList filelist = dir.entryInfoList();
    qDebug() << filelist;

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
    QString co1 = FindCommonString(candidateList);
    QString co2 = FindCommonString(candidateResList);
    qDebug() << " co1 = " << co1;
    qDebug() << " co2 = " << co2;

    if ( !co2.isEmpty())
    {
        this->CommandCompletion(co2);
    }
}
void Terminal::Analysis(const QString &cmdString)
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

QString Terminal::GetCommandString()
{
    QTextCursor cursor = ::textEdit->textCursor();
    QTextBlock textBlock = cursor.block();

    QString selectLine = textBlock.text();
    QString currStr = QDir::currentPath() + ">";
    QString cmdStr = selectLine.mid(currStr.length());
    QString trimedCmd = cmdStr.trimmed();
    qDebug() << "cmdStr="<<cmdStr;
    qDebug() << "trimedCmd="<<trimedCmd;
    return trimedCmd;
}

void Terminal::CommandCompletion(const QString &str)
{
    QTextCursor cursor = ::textEdit->textCursor();
    cursor.insertText( str );
    ::textEdit->setTextCursor(cursor);
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
    if ( strlist.length() == 0 ) return "";
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
