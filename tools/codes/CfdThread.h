#pragma once

#include <QObject>
#include <QThread>

class CfdThread : public QThread
{
    Q_OBJECT
public:
    CfdThread();
public:
    void beginCFD();
    void pauseCFD();
    void stopCFD();
private:
    void download( const std::string & file );
signals:
    void newValue(int seq, int diceValue);
protected:
    void run() override;
private:
    bool stop = false;
    bool paused = false;
};
