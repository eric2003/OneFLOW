/**
 * 
 */
function click_vis(){
	
	var t=document.getElementById("t").value;
	
	if(isNaN(t)||t==""){
		alert("温度值设置错误！！");
	}else{
		var tem=parseFloat(t);
		var miu0=parseFloat(document.getElementById("miu0").value);
		var t0=parseFloat(document.getElementById("t0").value);
		var ts=parseFloat(document.getElementById("ts").value);
		
		var miu=miu0*Math.pow(tem/t0,1.5)*(t0+ts)/(tem+ts);
		
		var mout=document.getElementById("miu-out");
		mout.value=miu.toExponential(4);
	}
}

function click_sound(){
	
	var t=document.getElementById("t").value;
	
	if(isNaN(t)||t==""){
		alert("温度值设置错误！！");
	}else{
		var tem=parseFloat(t);		
		var gama=1.4
		var R=287.2
		
		var c=Math.sqrt(gama*R*tem)
		document.getElementById("csound").value=c.toFixed(4);
		
		var mach=document.getElementById("machNum").value
		document.getElementById("speed").value=(mach*c).toFixed(4);
	}
}

function click_copyvis(){
		var mout=document.getElementById("miu-out");		
		document.getElementById("miu").value=mout.value;
}

function click_re(){
	
	var rho=parseFloat(document.getElementById("rho").value);
	var u=parseFloat(document.getElementById("u").value);
	var len=parseFloat(document.getElementById("len").value);
	var miu=parseFloat(document.getElementById("miu").value);
	
	var re=rho*u*len/miu;
	
	var reout=document.getElementById("re-out");
	reout.value=re.toExponential(2);
}