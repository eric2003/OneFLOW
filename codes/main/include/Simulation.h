/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
This file is part of OneFLOW.

OneFLOW is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OneFLOW is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#pragma once
#include "Configure.h"
#include <vector>

BeginNameSpace( ONEFLOW )

//调用仿真求解模块，进行求解和后处理的类
class Simulation
{
public:
	Simulation( int argc, char ** argv );
	virtual ~Simulation();
public:
	void Run();

public:
	void PreProcess();
	void MainProcess();
	void PostProcess();

protected:
	//初始化
	void InitSimu();
	//执行求解
	void RunSimu();
	
private:
	//命令行参数
	std::vector<std::string> args;
};

EndNameSpace
