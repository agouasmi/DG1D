#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <thread>

#include "discpp.h"

/*  ---  Functions to make Dislin more compact  --  */

void plot_init_setup(Dislin &g){
	
	/* (1) Setting of page format, file format and filename */

	g.scrmod("revers");
	g.sclmod("full");
	g.page(800,600);
	g.window(0,0,800,600);
	g.metafl("cons");

	/* (2) Initialization */

	g.disini();

};

void plot_axis_setup(Dislin &g){

	/* (3) Setting of plot parameters */ 

	g.hname(15);
	g.namdis(10,"xy");
	g.chaspc(-0.01);
	g.height(15);	
	g.htitle(15);
	g.vkytit(-20);
	
	/* (4) Setting axis and graph parameters*/

	g.ax2grf();
	g.frame(3);	
	g.ticks(3,"XY");
	g.ticpos("revers","XY"); 
	g.simplx();  //

	g.name("x","X");
	g.name("h(x)","Y");

	g.axspos(100,500);
	g.axslen(600,400);
		
};

void set_title(Dislin &g, string title_name){

	g.titlin("Position", 3);
	g.titlin(title_name.c_str(),4);

};

/* --- Pausing function for animations --- */

void pause(double seconds){
	std::this_thread::sleep_for(std::chrono::milliseconds(int(1000*seconds)));	
};
	
