const doc_pages = new Set([
	"BioSANS2020.analysis.numeric.sample_points.html",
	"BioSANS2020.analysis.numeric.transform_data.html",
	"BioSANS2020.analysis.plotting.plot_traj.html",
	"BioSANS2020.BioSANS.html",
	"BioSANS2020.BioSSL.html",
	"BioSANS2020.cli_functs.ssl_calls.html",
	"BioSANS2020.gui_functs.draw_figure.html",
	"BioSANS2020.gui_functs.prepare_canvas.html",
	"BioSANS2020.gui_functs.scrollable_text.html",
	"BioSANS2020.math_functs.sbml_math.html",
	"BioSANS2020.model.fileconvert.convtopotosbml.html",
	"BioSANS2020.model.fileconvert.process_sbml.html",
	"BioSANS2020.model.new_file.html",
	"BioSANS2020.model.ode_parse.ode_extract.html",
	"BioSANS2020.model.param_est.my_mcem.html",
	"BioSANS2020.model.param_est.param_estimate.html",
	"BioSANS2020.model.param_est.param_slider.html",
	"BioSANS2020.model.topology_view.html",
	"BioSANS2020.myglobal.mglobals.html",
	"BioSANS2020.myglobal.proc_global.html",
	"BioSANS2020.prepcodes.process.html",
	"BioSANS2020.prepcodes.processes_hub.html",
	"BioSANS2020.propagation.create_wxmaxima_command.html",
	"BioSANS2020.propagation.deterministic.euler_mod.html",
	"BioSANS2020.propagation.deterministic.lna_approx.html",
	"BioSANS2020.propagation.deterministic.lna_function_of_time.html",
	"BioSANS2020.propagation.deterministic.ode_int.html",
	"BioSANS2020.propagation.deterministic.runge_kutta4.html",
	"BioSANS2020.propagation.law_of_localization.html",
	"BioSANS2020.propagation.propensity.html",
	"BioSANS2020.propagation.recalculate_globals.html",
	"BioSANS2020.propagation.stochastic.gillespie_ssa.html",
	"BioSANS2020.propagation.stochastic.mystiffcle.html",
	"BioSANS2020.propagation.stochastic.mytauleap.html",
	"BioSANS2020.propagation.stochastic.tau_leaping.html",
	"BioSANS2020.propagation.stochastic.tau_leaping2.html",
	"BioSANS2020.propagation.symbolic.analytical_sol.html",
	"BioSANS2020.propagation.symbolic.lna_approx2.html",
	"BioSANS2020.RunBioSANS.html"				
]);

$.fn.extend({
	wru: function(){
		return $(this).each(function(){
			$(this).addClass("encapsulate");
		});
	}
});	
		
$(document).ready(function(){
	

	MyIFrame = $(".modal-body").find(".modal-iframe")

	$("li.toctree-l2").find("li.toctree-l2").find('a').siblings('ul').toggle();
	$("li.toctree-l2").find('a').click(function(event){
		event.preventDefault();
		$(this).siblings('ul').toggle();
	});

	$(".pre").click(function(event){
		event.preventDefault();
		file_html = $(this).text()+".html";
		if (doc_pages.has(file_html)){
			MyIFrame.attr("src", "pydoc_out/"+file_html)
		}else if(file_html == "Sphinx base documentation.html"){
			MyIFrame.attr("src", "https://efajiculay.github.io/biosans.docs/");
		}
	});			
				
	$(".window_installer").click(function(e) {
		e.preventDefault(); 
		file = $(this).text()
		if (file != "README.html"){
			window.location.href = "Installers/windows/"+file;
		}else{
			MyIFrame.attr("src", "Installers/windows/"+file);
		}
	});
	
	$(".ubuntu_installer").click(function(e) {
		e.preventDefault(); 
		file = $(this).text()
		if (file != "README.html"){
			window.location.href = "Installers/ubuntu/"+file;
		}else{
			MyIFrame.attr("src", "Installers/ubuntu/"+file);
		}
	});
	
	$(".macos_installer").click(function(e) {
		e.preventDefault(); 
		file = $(this).text()
		if (file != "README.html"){
			window.location.href = "Installers/macos/"+file;
		}else{
			MyIFrame.attr("src", "Installers/macos/"+file);
		}
	});
	
	$("#windows_id").toggle();
	$("#ubuntu_id").toggle();
	$("#macos_id").toggle();
	
	$("#git_hub").click(function(e){
		e.preventDefault(); 
		window.location.href = $(this).attr("href");
	});
	
	$("#test_pypi").click(function(e){
		e.preventDefault(); 
		window.location.href = $(this).attr("href");
	});
	
	$("#Stut_1").click(function(event){
		event.preventDefault();
		file_html = "https://www.youtube.com/embed/qt3e6Uw5F2g";
		MyIFrame.attr("src", file_html);
	});
	
	$("#Stut_2").click(function(event){
		event.preventDefault();
		file_html = "https://www.youtube.com/embed/nBODhuZ2aSg";
		MyIFrame.attr("src", file_html);
	});
	
	$("#Dtut_1").click(function(event){
		event.preventDefault();
		file_html = "https://www.youtube.com/embed/3DPYf5s6BRI";
		MyIFrame.attr("src", file_html);
	});
	
	$("#Ttut_1").click(function(event){
		event.preventDefault();
		file_html = "https://www.youtube.com/embed/Zxg7XfEAbPA";
		MyIFrame.attr("src", file_html);
	});
	
	$("#Ttut_2").click(function(event){
		event.preventDefault();
		file_html = "https://www.youtube.com/embed/A3JIYJZqHkI";
		MyIFrame.attr("src", file_html);
	});
	
				
	MyIFrame.on('load',function(){
		what = $(this).attr('src')
		if(what.indexOf("youtube") != -1){
			$(this).height('100vh');
		}else{
			$(this).height(10000);
		}
	});
											
	
});
