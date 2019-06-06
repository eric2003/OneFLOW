/**
 * 2018.08.15
 */

$(document).ready(function() {

	// 当算例为motion时，流动状态只能为transient
	setMotionFlow()

	$("#case_type").change(function() {
		setMotionFlow()
	})

	// order blend
	setOrderBlendDiv()

	$("#order_blend").change(function() {
		setOrderBlendDiv()
	})

	// flow_state
	setFlowStateDiv();

	$("#flow_state").change(function() {
		setFlowStateDiv();
	});

	// venka
	setVenkaDiv()

	$("#limiter").change(function() {
		setVenkaDiv()
	});
	
	//输出变量参数
	setVarDiv("FieldVars")
	$("#FieldVars").change(function() {
		setVarDiv("FieldVars")
	});
	
	setVarDiv("BoundVars")
	$("#BoundVars").change(function() {
		setVarDiv("BoundVars")
	});

	// 动网格参数
	setDyna()
	$("#dynaFlag").change(function() {
		setDyna()
	});
	
	setDynaMethod()
	$("#dyMeshMethod").change(function() {
		setDynaMethod()
	});

	// 添加运动组
	document.getElementById("NumMotionGroup").value=0
	$("#MotionGroupFlag").change(function() {
		//setMotionGroup()
	});

	// 折叠
	$("#fold").accordion();
})

// //////////////////////////////
function setDyna() {
	
	//检查动网格是否勾选
	var flag = document.getElementById("dynaFlag").checked;

	if (flag) {
		$("#dynamics").css("display", "block")
		document.getElementById("btnMgrp").disabled=false
		
		var ss = $("#flow_state").children('option:selected').val();
		if (ss == "steady") {
			document.getElementById("flow_state").value="transient";
			document.getElementById("case_type").value="motion";
			
			document.getElementById("time").value="dual-step";
				
			$("#steady").css("display", "none")
			$("#unsteady").css("display", "block")
		}		
	} else {
		$("#dynamics").css("display", "none")
		$("#MotionGroup").empty()
		document.getElementById("NumMotionGroup").value=0
		document.getElementById("btnMgrp").disabled=true
	}
}

function setDynaMethod(){
	//根据动网格方法显示不同的动网格参数
	var dym=document.getElementById("dyMeshMethod").value
	
	var dymethods=["#m_spring_origin","#m_rbf_origin",
		"#m_rbf_greedy","#m_rbf_idw"]
	
	dymethods.map(function(value) {
		$(value).css("display", "none")
	})
	
	if (dym=="spring-origin"||dym=="spring-wall") {
		$("#m_spring_origin").css("display", "block")
	}
	else if (dym=="rbf-origin") {
		$("#m_rbf_origin").css("display", "block")
	}
	else if (dym=="rbf-greedy") {
		$("#m_rbf_origin").css("display", "block")
		$("#m_rbf_greedy").css("display", "block")
	}
	else if (dym=="rbf-idw") {
		$("#m_rbf_origin").css("display", "block")
		$("#m_rbf_idw").css("display", "block")
	}
}

function resetMotionGroup() {
	$("#MotionGroup").empty()
	document.getElementById("NumMotionGroup").value=0
}

// motionGroup
function StrMotionType(n) {
	var motiontype = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">运动类型</span> <select id=\"motiontype"
			+ n + "\"" + "class=\"form-control\">" 
			+ "<option>dof3</option>"
			+ "<option>dof6</option>" 
			+ "<option>attitude</option>"
			+ "<option>pitch</option>" 
			+ "<option>udf-compiled</option>"
			+ "<option>udf-interpreted</option>"
			+ "<option>udf-deform</option>" + "</select>" + "</div>"

	return motiontype
}
function StrUDF(n) {
	var udf = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">UDF文件名称</span><input type=\"text\""
			+ "value=\"null\" class=\"form-control\" id=\"udfFileName" + n
			+ "\">"
			+ "<span class=\"input-group-addon\">UDF函数名称</span><input type=\"text\""
			+ "value=\"null\" class=\"form-control\" id=\"udfFuncName" + n
			+ "\">" + "</div>"

	return udf
}
function StrMass(n) {
	var udf = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">质量</span><input type=\"text\""
			+ "value=\"0\" class=\"form-control\" id=\"Mass" + n + "\">"
			+ "</div>"

	return udf
}
function StrVec(val, n) {
	var str = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">" + val + "</span>"
			+ "<span class=\"input-group-addon\">x</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-x" + "\">"
			+ "<span class=\"input-group-addon\">y</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-y" + "\">"
			+ "<span class=\"input-group-addon\">z</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-z" + "\">" + "</div>"

	return str
}

function StrMat(val, n) {
	var str = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">" + val + "</span>"
			+ "<span class=\"input-group-addon\">xx</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-xx" + "\">"
			+ "<span class=\"input-group-addon\">yy</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-yy" + "\">"
			+ "<span class=\"input-group-addon\">zz</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-zz" + "\">"
			
			+ "<span class=\"input-group-addon\">xy</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-xy" + "\">"
			+ "<span class=\"input-group-addon\">xz</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-xz" + "\">"
			+ "<span class=\"input-group-addon\">yz</span>"
			+ "<input value=\"0\" type=\"text\" class=\"form-control\" id=\""
			+ val + n + "-yz" + "\">" + "</div>"

	return str
}

function StrHull(n) {
	var udf = "<div class=\"input-group\">"
			+ "<span class=\"input-group-addon\">是否是洞点<input type=\"checkbox\""
			+ "id=\"isHullPoint" + n + "\"></span>" + "</div>"

	return udf
}
function setMotionGroup() {

	// 先清空
	//$("#MotionGroup").empty()

	//var flag = document.getElementById("MotionGroupFlag").checked;
	var dflag = document.getElementById("dynaFlag").checked;
	var i = document.getElementById("NumMotionGroup").value;
	++i;

	if (dflag) {
			var str = "<div id=\"Object-" + i + "\">" + "Object-" + i
					+ "</div>";
			
			$("#MotionGroup").append(str)
			$("#MotionGroup").append(StrMotionType(i))
			$("#MotionGroup").append(StrUDF(i))
			$("#MotionGroup").append(StrMass(i))
			$("#MotionGroup").append(StrMat("惯性矩阵", i))
			$("#MotionGroup").append(StrVec("质心", i))
			$("#MotionGroup").append(StrVec("欧拉角", i))
			$("#MotionGroup").append(StrVec("线速度", i))
			$("#MotionGroup").append(StrVec("角速度", i))
			$("#MotionGroup").append(StrHull(i))
	}
	
	document.getElementById("NumMotionGroup").value=i;
}

function setMotionFlow() {
	var ct = document.getElementById("case_type").value

	if (ct == "motion") {
		document.getElementById("flow_state").value = "transient"
		setFlowStateDiv()
	}
}

function setFlowStateDiv() {
	var ss = $("#flow_state").children('option:selected').val();
	if (ss == "steady") {
		$("#steady").css("display", "block")
		$("#unsteady").css("display", "none")
	} else {
		$("#steady").css("display", "none")
		$("#unsteady").css("display", "block")
	}
}

function setOrderBlendDiv() {
	if (document.getElementById("order_blend").checked) {
		$("#ord_blend").css("display", "block")
	} else {
		$("#ord_blend").css("display", "none")
	}
}

function setVenkaDiv() {
	var ss = $("#limiter").children('option:selected').val()

	if (ss == "Venkatakrishnan") {
		$("#venka").css("display", "block")
	} else {
		$("#venka").css("display", "none")
	}
}

//设置流量变量可见性
function setVarDiv(name) {
	var flag = document.getElementById(name).checked;

	if (flag) {
		$("#"+name+"List").css("display", "block")
	} else {
		$("#"+name+"List").css("display", "none")
	}
}

function setHide(name) {
	$("#"+name).css("display", "none")
}
function setShow(name) {
	$("#"+name).css("display", "block")
}

// 复制到剪切板
function copyParam() {
	var obj = document.getElementById("jsonContent");
	obj.select(); // 选择对象
	document.execCommand("Copy"); // 执行浏览器复制命令
	alert("已复制");
}

function generate() {

	var cjson = {

	}

	// 元数据
	var jmeta = {
		"CodeName" : "Caldron",
		"InputFmtVersion" : "2.6",
		"KernelVersion" : "2.0"
	}

	// 基本参数
	var jbasic = {
		"case.name" : document.getElementById("case.name").value,
		"case.type" : document.getElementById("case_type").value,
		"flow.state" : document.getElementById("flow_state").value,
		"case.method" : document.getElementById("case.method").value
	}

	// 格式和相关参数
	var jscheme = {
		"flux" : document.getElementById("flux").value,
		"time" : document.getElementById("time").value,
		"gradient" : document.getElementById("gradient").value,
		"limiter" : document.getElementById("limiter").value
	}

	var lim = $("#limiter").children('option:selected').val()

	if (lim == "Venkatakrishnan") {
		jscheme["Venkatakrishnan-K"] = parseFloat(document
				.getElementById("venkaK").value)
	}
	
	jscheme["limiterStencil"] = document.getElementById("limiterStencil").value

	// 参考量
	var x = document.getElementById("refPx").value
	var y = document.getElementById("refPy").value
	var z = document.getElementById("refPz").value

	var jreference = {
		"ref.length" : parseFloat(document.getElementById("refLength").value),
		"ref.area" : parseFloat(document.getElementById("refArea").value),
		"ref.point" : [ parseFloat(x), parseFloat(y), parseFloat(z) ]
	}

	// 求解控制
	var jcontrol = {
		"step.begin" : parseInt(document.getElementById("step.begin").value),
		"step.end" : parseInt(document.getElementById("step.end").value),
		"step.put" : parseInt(document.getElementById("step.put").value),
		"initByMapping" : document.getElementById("initByMapping").checked,
		"visciousity" : document.getElementById("visciousity").value
	}

	var order = document.getElementById("space.order").value
	if (order == "first order") {
		jcontrol["space.order"] = 1
	} else {
		jcontrol["space.order"] = 2
	}

	if (document.getElementById("order_blend").checked) {
		var blend0 = document.getElementById("blend0").value
		var blend1 = document.getElementById("blend1").value
		jcontrol["order.blend"] = [ parseInt(blend0), parseInt(blend1) ]
	} else {
		jcontrol["order.blend"] = [ 0, 0 ]
	}

	jcontrol["EnableLocalStep"] = document.getElementById("EnableLocalStep").checked
	jcontrol["EnableResSmooth"] = document.getElementById("EnableResSmooth").checked

	if (jbasic["flow.state"] == "steady") {
		var cfl0 = document.getElementById("cfl_num0").value
		var cfl1 = document.getElementById("cfl_num1").value

		var cflramp0 = document.getElementById("cfl_ramp0").value
		var cflramp1 = document.getElementById("cfl_ramp1").value

		var jsteady = {
			"res.end" : parseFloat(document.getElementById("res_end").value),
			"cfl.num" : [ parseFloat(cfl0), parseFloat(cfl1) ],
			"cfl.ramp" : [ parseInt(cflramp0), parseInt(cflramp1) ]
		}

		jcontrol["Steady"] = jsteady
	} else {

		// 时间步长
		var jstep = {
			"EnableFixedStep" : document.getElementById("usefixstep").checked,
			"FixedStepSize" : parseFloat(document.getElementById("step_size").value)
		};

		// 双时间步长参数
		var jDts = {
			"Scheme" : "LU-SGS",
			"MaxIter" : parseInt(document.getElementById("inner_iter").value),
			"ResEnd" : parseFloat(document.getElementById("inner_res").value),
			"CflNum" : parseFloat(document.getElementById("inner_cfl").value)
		};

		// 非定常参数
		var jtransient = {
			"TimeStep" : jstep,
			"DualTimeStep" : jDts
		}

		jcontrol["Transient"] = jtransient
	}

	// out-ctrl
	var jout = {
		"flowfilefmt" : document.getElementById("flowfilefmt").value, // 流场文件格式{vtk,tecplot}
		"forcefilefmt" : "csv", // 气动力文件格式{csv,tecplot}
		"residualfilefmt" : "csv", // 残差文件格式{csv,tecplot}
		"enableSpartanTecplot" : false
	}
	
	//要输出的流场变量
	if(document.getElementById("FieldVars").checked){	
		jout["FieldVar"]=[]
		for (var i = 0; i < 13; i++) {
			var e=document.getElementById("fv"+i);
			if (e.checked) {
				jout["FieldVar"].push(e.value)
			}
		}
	}
	else {
		jout["FieldVar"]=[ "ro", "u", "v", "w", "p", "T", "Ma" ]; // 要输出的流场变量
	}
	
	// 要输出的边界变量
	if(document.getElementById("BoundVars").checked){
		jout["BoundVar"]=[]
		for (var i = 0; i < 15; i++) {
			var e=document.getElementById("bv"+i);
			if (e.checked) {
				jout["BoundVar"].push(e.value)
			}
		}
	}
	else {
		jout["BoundVar"]=[ "ro", "u", "v", "w", "p", "T", "Ma" ]; // 要输出的流场变量
	}

	cjson["meta"] = jmeta
	cjson["basic"] = jbasic
	cjson["control"] = jcontrol
	cjson["scheme"] = jscheme
	cjson["reference"] = jreference
	cjson["out-ctrl"] = jout

	// 动网格参数
	var dyna = document.getElementById("dynaFlag").checked;

	if (dyna) {
		jdyna = {
			"Method" : document.getElementById("dyMeshMethod").value
		}
		
		//网格重构
		jdyna["EnableRemesh"]=document.getElementById("EnableRemesh").checked;
		//网格质量检测标准
		jdyna["MeshQualityStandard"]=document.getElementById("MeshQualityStandard").value;
		
		//弹簧方法
		if (jdyna["Method"]=="spring-origin"||jdyna["Method"]=="spring-wall") {
			jdyna["Spring.MaxIteration"]=parseInt(
					document.getElementById("springIter").value);
			jdyna["Spring.EnableTorsion"]=document.getElementById("springTorsion").checked
		}
		
		//RBF类方法
		if (jdyna["Method"]=="rbf-origin") {
			jdyna["rbf.SupportRadius"]=parseFloat(document.getElementById("rbfSupportRadius").value);
		}
		
		if (jdyna["Method"]=="rbf-greedy"||jdyna["Method"]=="rbf-idw") {
			jdyna["rbf.SupportRadius"]=parseFloat(document.getElementById("rbfSupportRadius").value);
			jdyna["rbf.GreedyError"]=parseFloat(document.getElementById("rbfGreedyError").value);
			jdyna["rbf.GreedyinitPts"]=parseInt(document.getElementById("rbfGreedyinitPts").value);
		}
		
		//聚类RBF方法
		if (jdyna["Method"]=="rbf-idw") {
			jdyna["rbf.ClusterInitK"]=parseInt(document.getElementById("rbfClusterInitK").value);	
		}
		
		cjson["dynamics"] = jdyna
		
		//运动组
		jobjects={};
		
		var nObj=document.getElementById("NumMotionGroup").value

		for (var i = 1; i <= nObj; i++) {
			var jobj="Object"+i;
			
			//motion type
			jobjects[jobj]={
					"Motion.type":document.getElementById("motiontype"+i).value,
			}
			
			//udf
			var mt=document.getElementById("motiontype"+i).value
			if (mt=="udf-compiled"||mt=="udf-interpreted"||mt=="udf-deform") {
				jobjects[jobj]["UDF"]={
						"File.name":document.getElementById("udfFileName"+i).value,
						"Function.name":document.getElementById("udfFuncName"+i).value
				}
			}
			
			//motion para
			var jmpara={}
			jmpara["Mass"]=parseFloat(document.getElementById("Mass"+i).value)

			var xx=parseFloat(document.getElementById("惯性矩阵"+i+"-xx").value)
			var yy=parseFloat(document.getElementById("惯性矩阵"+i+"-yy").value)
			var zz=parseFloat(document.getElementById("惯性矩阵"+i+"-zz").value)
			var xy=parseFloat(document.getElementById("惯性矩阵"+i+"-xy").value)
			var xz=parseFloat(document.getElementById("惯性矩阵"+i+"-xz").value)
			var yz=parseFloat(document.getElementById("惯性矩阵"+i+"-yz").value)
			
			jmpara["InertialMatrix"]=[
				[xx,xy,xz],
				[xy,yy,yz],
				[xz,yz,zz]
			]
			
			var cenx=parseFloat(document.getElementById("质心"+i+"-x").value)
			var ceny=parseFloat(document.getElementById("质心"+i+"-y").value)
			var cenz=parseFloat(document.getElementById("质心"+i+"-z").value)
			
			jmpara["MassCenter"]=[cenx,ceny,cenz]
			
			//欧拉角
			var ax=parseFloat(document.getElementById("欧拉角"+i+"-x").value)
			var ay=parseFloat(document.getElementById("欧拉角"+i+"-y").value)
			var az=parseFloat(document.getElementById("欧拉角"+i+"-z").value)
			
			jmpara["EulerAngle"]=[ax,ay,az]
			
			//线速度
			var vx=parseFloat(document.getElementById("线速度"+i+"-x").value)
			var vy=parseFloat(document.getElementById("线速度"+i+"-y").value)
			var vz=parseFloat(document.getElementById("线速度"+i+"-z").value)
			
			jmpara["LinearVelocity"]=[vx,vy,vz]
			
			//角速度
			var wx=parseFloat(document.getElementById("角速度"+i+"-x").value)
			var wy=parseFloat(document.getElementById("角速度"+i+"-y").value)
			var wz=parseFloat(document.getElementById("角速度"+i+"-z").value)
			
			jmpara["AngularVelocity"]=[wx,wy,wz]
			
			jobjects[jobj]["MotionPara"]=jmpara
			
			//hull point
			jobjects[jobj]["isVertexHull"]=document.getElementById("isHullPoint"+i).checked
		}
		
		cjson["MotionGroup"]=jobjects;
	}

	var x = document.getElementById("jsonContent")
	x.value = JSON.stringify(cjson, null, 4);
}
