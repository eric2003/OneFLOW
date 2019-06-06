
//设置网格参数


var nsurfs=0;	//边界的数目

$(document).ready(function() {
	resetPartGroup();	
})

function setTw()
{
	//等温壁面
	for (var i = 0; i < nsurfs; i++) {
		var attr=document.getElementById("surfAttr"+i).value

		alert(attr)

		$("#surfAttr"+i).change(function() {
			if(attr=="wall-isotherm"){
				$("#defaultTw").empty()
				var defTw=0 //jm["BoundaryConditions"]["default.WallTemperature"]
				var str=Span("默认值","默认值")+Input("defTw",defTw)
				$("#defaultTw").append(str);				
			}
			//break;
		});
	}
}

//复制到剪切板
function copyMeshParam() {
	var obj = document.getElementById("meshJson");
	obj.select(); // 选择对象
	document.execCommand("Copy"); // 执行浏览器复制命令
	alert("已复制");
}

function Input(id,value){
	return "<input id=\""+id+"\" value=\""
	+value+
	"\" class=\"form-control\" type=\"text\">";
}
function Span(id,value){
	return "<span id=\""+id+"\""+" class=\"input-group-addon\">"+value+"</span>"
}
function HiddenInput(id,value){
	return "<input type=\"hidden\" id=\""+id+"\" value=\""
	+value+
	"\" class=\"form-control\" type=\"text\">";
}
function resetPartGroup() {
	$("#partGroup").empty()
	$("#partAll").empty()
	document.getElementById("NumPartGroup").value=0
}
function setPartGroupValue(i) {
	var jmesh = document.getElementById("meshJson").value;		
	var jm=JSON.parse(jmesh);//["Mesh"]

	//surface
	var surfs=jm["Mesh"]["Surface"]
	s="";
	var count=0;
	for (var name in surfs) {	
		if (name.substring(0,2)=="__") {
			continue;
		}

		var id="pgrp_"+name;
		var chked=document.getElementById(id).checked

		if (chked) {
			if (count>0) {
				s+=","+name;
			}
			else {
				s+=name;
			}
			++count;
			//记录之后清除勾选
			document.getElementById(id).checked=false;
		}
	}

	if (count==0) {
		alert("请勾选下面的边界名称")
	}
	else {
		document.getElementById("pgrpList_"+i).value=s;
	}	
}
function setPartGroup() {

	var i =document.getElementById("NumPartGroup").value;

	if (i==0) {
		var jmesh = document.getElementById("meshJson").value;		
		var jm=JSON.parse(jmesh);//["Mesh"]

		//surface
		var surfs=jm["Mesh"]["Surface"]
		for (var name in surfs) {	
			if (name.substring(0,2)=="__") {
				continue;
			}
			var s="<div  class=\"input-group\">"+
			"<input type=\"checkbox\" id=\"pgrp_"
			+name+"\">"+
			"<span class=\"input-group-addon\">"
			+name+"</span></div>"
			$("#partAll").append(s)
		}
	}

	++i;
	str1="<div  class=\"input-group\"><span class=\"input-group-addon\">"+
	"<button id=\"btnPgrp\" type=\"button\" class=\"btn btn-primary\" onclick=\"setPartGroupValue("+i+")\">设置</button>"+
	"</span>";

	var str = "<input type=\"text\" id=\"pgrp_"+i+"\" value=\"Group"+
	i+
	"\" class=\"form-control\" type=\"text\">"+
	"<input class=\"form-control\" type=\"text\" id=\"pgrpList_"+i+"\">"+
	"</div>";

	$("#partGroup").append(str1+str)

	document.getElementById("NumPartGroup").value=i;
}
function initMeshPara() {

	var jmesh = document.getElementById("meshJson").value;

	if(jmesh==null||jmesh==""){
		alert("输入参数不能为空")
	}

	var jm=JSON.parse(jmesh);//["Mesh"]

	var scale=jm["Mesh"]["Scale"]

	$("#surfs").empty()
	$("#bodies").empty()
	resetPartGroup()

	//缩放
	document.getElementById("scaleX").value=scale[0]
	document.getElementById("scaleY").value=scale[1]
	document.getElementById("scaleZ").value=scale[2]

	var str1="<div class=\"input-group\">"
		var strspan="<span class=\"input-group-addon\">"

			//body
			var bodies=jm["Mesh"]["Body"]

	var i=0	
	for (var name in bodies) {

		var media="<select id=\"bodyMedia"+i+"\""+
		"class=\"form-control\">"+
		"<option>fluid</option>"+
		"<option>solid</option>"+
		"</select>"

		var part=bodies[name]["Id"]

		var str=str1+
		Span("bodyName"+i,name)+media+
		HiddenInput("bodyPart"+i,part)+
		"</div>"

		$("#bodies").append(str)
		$("#boco"+i).val(bodies[name]["Media"])
		++i
	}
	$("#bodies").append(Span("nbodies",i))

	//surface
	var surfs=jm["Mesh"]["Surface"]

	i=0;	
	var countHidden=0;
	var countIsothermWall=0;

	for (var name in surfs) {		
		var boco="<select id=\"surfAttr"+i+"\""+
		"class=\"form-control\">"+
		"<option>inflow</option>"+
		"<option>outflow</option>"+
		"<option>farfield</option>"+
		"<option>symmetry</option>"+
		"<option>wall-adiabatic</option>"+
		"<option>wall-isotherm</option>"+
		"<option>wall-slip</option>"+
		"<option>none</option>"+
		"</select>"

		var mgrp=surfs[name]["MotionGroup"]

		if (mgrp==undefined) {
			mgrp=-1
		}

		var part=surfs[name]["Id"]		
		var dom=surfs[name]["Domain"]		

		var str="";
		var strHid="<div class=\"input-group\" id=\""+"surfDiv_"+name+"\">";

		if(dom==null||dom==undefined){
			str=strHid+
			Span("surfName"+i,name)+boco
			+HiddenInput("surfPart"+i,part)
			+Span("span_motionGroup"+i,"MotionGroup")
			+Input("MotionGroup"+i,mgrp)
			+"</div>";

			$("#surfs").append(str);
			$("#surfAttr"+i).val(surfs[name]["Attr"]);
		}
		else {
			str=strHid
			+Span("surfName"+i,name)+boco
			+HiddenInput("surfPart"+i,part)
			+HiddenInput("surfDomainL"+i,dom[0])
			+HiddenInput("surfDomainR"+i,dom[1])
			+Span("span_motionGroup"+i,"MotionGroup")
			+Input("MotionGroup"+i,mgrp)
			+"</div>";

			$("#surfs").append(str);
			$("#surfDiv_"+name).hide();
			$("#surfAttr"+i).val(surfs[name]["Attr"]);
			++countHidden;
		}

		//等温壁面
		var bc=surfs[name]["Attr"]
		if(bc=="wall-isotherm"){
			++countIsothermWall;
		}

		++i
	}

	nsurfs=i;
	$("#surfs").append(Span("nsurfs",i))

	//默认来流值
	$("#defaultInflow").empty()
	var defInflow=jm["BoundaryConditions"]["default.InflowParameter"]

	var str=Span("默认值","默认值")+Span("p","p")+Input("defp",defInflow["p"])+
	Span("T","T")+Input("defT",defInflow["t"])+
	Span("u","u")+Input("defu",defInflow["u"])+
	Span("v","v")+Input("defv",defInflow["v"])+
	Span("w","w")+Input("defw",defInflow["w"])

	$("#defaultInflow").append(str)

	//默认壁温
	if(countIsothermWall>0){
		$("#defaultTw").empty()
		var defTw=jm["BoundaryConditions"]["default.WallTemperature"]
		var str=Span("默认值","默认值")+Input("defTw",defTw)
		$("#defaultTw").append(str)
	}		
}

function generateMeshJson(){

	//
	var nsurf=document.getElementById("nsurfs").innerText
	var jsurf={}

	var nIsothermWall=0;
	for (var i = 0; i < nsurf; i++) {
		var name=document.getElementById("surfName"+i).innerText

		jsurf[name]={
			"Attr":document.getElementById("surfAttr"+i).value,
			"Id":parseInt(document.getElementById("surfPart"+i).value),
			"MotionGroup":parseInt(document.getElementById("MotionGroup"+i).value)
		}

		var domL=document.getElementById("surfDomainL"+i)
		var domR=document.getElementById("surfDomainR"+i)
		if (domL!=undefined&&domL!=null) {
			jsurf[name]["Domain"]=[parseInt(domL.value),parseInt(domR.value)]
		}

		//等温壁面
		if(jsurf[name]["Attr"]=="wall-isotherm"){
			++nIsothermWall;
		}
	}

	var nbody=document.getElementById("nbodies").innerText
	var jbody={}

	for (var i = 0; i < nbody; i++) {
		jbody[document.getElementById("bodyName"+i).innerText]={
				"Id":parseInt(document.getElementById("bodyPart"+i).value),
				"Media":document.getElementById("bodyMedia"+i).value			
		}
	}
	//缩放
	var x = document.getElementById("scaleX").value
	var y = document.getElementById("scaleY").value
	var z = document.getElementById("scaleZ").value

	//partGroup
	var jpgroup={}

	var npar=document.getElementById("NumPartGroup").value

	for (var i = 1; i <= npar; i++) {
		var jpid="pgrp_"+i;
		var name=document.getElementById(jpid).value;
		var id2="pgrpList_"+i
		var s=document.getElementById(id2).value;

		if (s.length!=0) {
			var a = s.split(',');
			jpgroup[name]=a;
		}		
	}

	//默认值
	jdef={
			"p":parseFloat(document.getElementById("defp").value),
			"t":parseFloat(document.getElementById("defT").value),
			"u":parseFloat(document.getElementById("defu").value),
			"v":parseFloat(document.getElementById("defv").value),
			"w":parseFloat(document.getElementById("defw").value)
	}
	
	var cjson = {}
	
	//边界条件
	cjson["BoundaryConditions"]={
			"default.InflowParameter":jdef
	};

	if (nIsothermWall>0) {
		var defTw=parseFloat(document.getElementById("defTw").value);
		cjson["BoundaryConditions"]["default.WallTemperature"]=defTw;
	}
	
	//网格信息
	jmesh={
			"Body":jbody,
			"Scale":[ parseFloat(x), parseFloat(y), parseFloat(z) ],
			"PartGroup":jpgroup,
			"Surface":jsurf
	}

	cjson["Mesh"]=jmesh

	document.getElementById("meshJson").value = JSON.stringify(cjson, null, 4);
}