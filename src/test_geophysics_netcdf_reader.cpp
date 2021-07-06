/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "mpi.h"
#include "netcdf.h"

#include <cstdio>
#include <netcdf>
#include <vector>
#include <limits>

#ifdef USEGLOBALSTACKTRACE
	#include "stacktrace.h"
	cStackTrace gtrace;
#endif

#include "marray.hxx"
using namespace andres;

using namespace netCDF;
using namespace netCDF::exceptions;

#include "general_utils.h"
#include "file_utils.h"
#include "vector_utils.h"
#include "file_formats.h"
#include "geophysics_netcdf.h"
#include "gdal_utils.h"
#include "stopwatch.h"
#include "logger.h"


class cLogger glog;

bool example_magnetics(){

	bool status;		
	std::string indir,ncpath;
	
	indir   = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/";
	ncpath  = indir + "GSSA_P1255MAG_Warrina.nc";
             
	//indir  = R"(Z:\projects\geophysics_netcdf\ncfiles\)";	
	//ncpath = indir + "P1152MAG_V2.nc";


	//Open the file and initialise the indexes
	cGeophysicsNcFile ncfile(ncpath, NcFile::read);

	//Get the line numbers using the standard_name attributes
	std::vector<int> linenumber   = ncfile.getLineNumbers();
	std::vector<int> flightnumber = ncfile.getFlightNumbers();

	std::string stdname_x         = "longitude";
	std::string stdname_y         = "latitude";
	std::string stdname_tielev    = "total_magnetic_intensity_anomaly_tie_line_levelled";
	std::string stdname_mlev      = "total_magnetic_intensity_anomaly_micro_levelled";
	std::string stdname_awagslev  = "total_magnetic_intensity_anomaly_datum_levelled";
	
	std::string xvarname = ncfile.getVarNameByLongName(stdname_x);
	std::string yvarname = ncfile.getVarNameByLongName(stdname_y);
	std::string mvarname = ncfile.getVarNameByLongName(stdname_awagslev);

	//Get 21th line number and index
	size_t lnum = linenumber[20];
	size_t lind = ncfile.getLineIndex(lnum);
	
	//Get mag data for 21th line by number and index
	std::vector<double> v1darray;
	status = ncfile.getDataByLineNumber(mvarname, lnum, v1darray);
	status = ncfile.getDataByLineIndex(mvarname, lind, v1darray);

	//Loop over all lines getting the x, y, and mag values
	std::vector<double> x, y; std::vector<float>  mag;
	cSampleVar xvar = ncfile.getSampleVar(xvarname);
	cSampleVar yvar = ncfile.getSampleVar(yvarname);
	cSampleVar mvar = ncfile.getSampleVar(mvarname);
	
	double t, ta, tb;
	ta = gettime(); status = xvar.getAll(x); tb = gettime();
	glog.logmsg("Get all of x (%lu doubles) - time=%lf s\n", x.size(), tb - ta);
	ta = gettime(); status = yvar.getAll(x); tb = gettime();
	glog.logmsg("Get all of y (%lu doubles) - time=%lf s\n", x.size(), tb - ta);
	ta = gettime(); status = mvar.getAll(mag); tb = gettime();
	glog.logmsg("Get all of mag (%lu floats) - time=%lf s\n", x.size(), tb - ta);

	t = 0;
	for (size_t li = 0; li < ncfile.nlines(); li += 100){
		ta = gettime();	
		status = xvar.getLine(li, x);
		status = yvar.getLine(li, y);
		status = mvar.getLine(li, mag);		
		tb = gettime();
		glog.logmsg("%lu nsamples=%lu - time=%lf s\n", li, x.size(), tb - ta);
		t += (tb - ta);
	}		
	glog.logmsg("Total every 100th line - time=%lf\n", t);
	return true;
}

bool example_aem_conductivity(){
	bool status; double t1, t2;
	//std::string indir  = R"(Z:\projects\geophysics_netcdf\ncfiles\)";	
	std::string indir = R"(Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\aem\ncfiles\)";

	//std::string indir  = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AEM_examples/";
	std::string ncpath = indir + "AUS_10008_WestK_LCI.nc";	
	//Open the file, get the line numbers
	cGeophysicsNcFile ncfile(ncpath, NcFile::write);						
	std::vector<int> linenumber = ncfile.getLineNumbers();
	//Determine conductivity variable name by its standard name attribute
	std::string stdname_conductivity = "layer_conductivity";
	std::string varname = ncfile.getVarNameByLongName(stdname_conductivity);	
	//Get conductivity variable
	cSampleVar vc = ncfile.getSampleVar(varname);
	//Get conductivity data all at once in 1d array
	std::vector<float> c1darray;
	t1 = gettime();	status = vc.getAll(c1darray); t2 = gettime();
	glog.logmsg("Get all at once (%lu floats) - time=%lf s\n", c1darray.size(), t2-t1);
	//Get conductivity line by line and band by band in 1d array
	t1 = gettime();
	for (size_t li = 0; li < ncfile.nlines(); li++){	
		for (size_t bi = 0; bi < vc.nbands(); bi++){			
			//status = vc.getLine(li, bi, c1darray);
			if (status == false)glog.logmsg("Error");
		}
	} t2 = gettime(); glog.logmsg("Get line by line and band by band in 1d array - time=%lf s\n", t2 - t1);
	//Get conductivity line by line in 2d array
	t1 = gettime();
	std::vector<std::vector<float>> c2darray;
	for (size_t li = 0; li < ncfile.nlines(); li++){	
		status = ncfile.getDataByLineIndex(varname, li, c2darray);
		if (status == false)glog.logmsg("Error");
	} t2 = gettime(); glog.logmsg("Get line by line in 2D array - time=%lf s\n", t2-t1);
	return true;	
}

bool test_create(){
	std::string ncpath = "test.nc";
	deletefile(ncpath);
	std::vector<size_t> linenumbers   = { 100, 200, 300, 400 };
	std::vector<size_t> linensamples  = { 10,   20,  30,  40 };
	std::vector<size_t> flightnumbers = {11, 22, 33, 44 };

	cGeophysicsNcFile   nc(ncpath,NcFile::FileMode::replace);
	nc.InitialiseNew(linenumbers, linensamples);
	size_t nwindows = 45;
	size_t nlayers  = 30;
	size_t nrxcomponents = 3;

	size_t ntotalsamples = nc.ntotalsamples();
	std::vector<int> fid = increment(ntotalsamples,0,1);
	std::vector<int> layers = increment(nlayers,0,1);
	std::vector<int> windows = increment(nwindows,0,1);
	std::vector<int> rxcomponents = increment(nrxcomponents,0,1);
	
	NcDim dim_rxcomponent = nc.addDimVar("rxcomponents", rxcomponents);
	NcDim dim_window = nc.addDimVar("windows", windows);	
	NcDim dim_layer  = nc.addDimVar("layers", layers);
	
	cSampleVar vfid = nc.addSampleVar("fiducial", ncInt);
	vfid.add_long_name("fiducial");
	vfid.add_units("1");
	vfid.add_missing_value(34);

	cSampleVar vx   = nc.addSampleVar("easting",  ncDouble);
	vx.add_long_name("X");
	vx.add_units("m");
	

	cSampleVar vy   = nc.addSampleVar("northing", ncDouble);	
	cSampleVar vconductivity = nc.addSampleVar("conductivity", ncDouble, dim_layer);
	cSampleVar vthickness    = nc.addSampleVar("thickness", ncDouble, dim_layer);
		
	vconductivity.add_long_name("conductivity");
	vconductivity.add_units("mS/m");
	vconductivity.add_missing_value(-999);
	
	vthickness.add_long_name("thickness");
	vthickness.add_units("m");
	vthickness.add_missing_value(-999);

	cLineVar vflight = nc.addLineVar("flight", ncInt);
	vflight.putAll(flightnumbers);

	std::vector<NcDim> emdims = { dim_rxcomponent, dim_window };
	cSampleVar vem = nc.addSampleVar("em", ncDouble, emdims);
		
	const size_t n = ntotalsamples*nrxcomponents*nwindows;
	std::vector<double> em = increment(n,0.0,1.0);		
	vfid.putAll(fid);
	vem.putAll(em);

	for (size_t li = 0; li < nc.nlines(); li++){
		size_t nls = nc.nlinesamples(li);
		std::vector<double> x = increment(nls,500000.0,10.0);
		std::vector<double> y = increment(nls,6500000.0,10.0);
		
		vx.putLine(li, x);
		vy.putLine(li, y);
		
		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> c(nls, li*10.0+bi);
			vconductivity.putLineBand(li, bi, c);
		}

		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> t(nls, li * 100.0 + bi);
			vthickness.putLineBand(li, bi, t);
		}
		
	}		
	return true;
};

bool test_update(){
	std::string ncpath = "test.nc";
	cGeophysicsNcFile   nc(ncpath, NcFile::FileMode::write);	
	size_t ntotalsamples = nc.ntotalsamples();
	cLineVar   lv = nc.addLineVar("extralinevar", ncDouble);
	cSampleVar sv = nc.addSampleVar("extrasamplevar", ncDouble);
	cSampleVar svw = nc.addSampleVar("extrasamplevarwindow", ncDouble,nc.getDim("windows"));

	NcDim ed = nc.addDim("extradim", 400);
	cSampleVar sve = nc.addSampleVar("extrawindowsamplevarextradim", ncDouble, ed);

	cSampleVar ev = nc.getSampleVar("easting");
	std::vector<double> e;
	bool status = ev.getAll(e);
	status = ev.getLine(1, e);
	e += 1.0;	
	status = ev.putLine(1, e);

	cSampleVar vem = nc.getSampleVar("em");
	std::vector<double> em;
	status = vem.getAll(em);
	status = vem.getLine(0,em);
	em *= 0.0;	
	status = vem.putLine(0,em);	
	return true;
};

bool test_aseggdfexport_1d(){
	//std::string indir  = R"(Z:\projects\geophysics_netcdf\ncfiles\)";
	std::string indir = R"(Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\awags_conversions\ncfiles\)";
	std::string ncpath  = indir + "P1152RAD.nc";
	std::string datpath = indir + "P1152RAD.dat";
	std::string dfnpath = indir + "P1152RAD.dfn";
	cGeophysicsNcFile nc(ncpath, NcFile::FileMode::read);
	nc.export_ASEGGDF2(datpath,dfnpath);
	return true;
};

bool test_aseggdfexport_2d(){
	//std::string indir = R"(Z:\projects\geophysics_netcdf\aem\)";
	std::string indir   = R"(Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\aem\ncfiles\)";
	std::string ncpath  = indir + "AUS_10008_WestK_LCI.nc";
	std::string datpath = indir + "AUS_10008_WestK_LCI.dat";
	std::string dfnpath = indir + "AUS_10008_WestK_LCI.dfn";
	cGeophysicsNcFile nc(ncpath, NcFile::FileMode::read);
	nc.export_ASEGGDF2(datpath, dfnpath);
	return true;
};

bool test_columnfile(){
	std::string datpath = R"(z:\projects\earth_sci_test\test_data\output\inversion.output.dat)";
	std::string dfnpath = R"(z:\projects\earth_sci_test\test_data\output\inversion.output.dfn)";

	cColumnFile A(datpath,dfnpath);
	bool status = A.readnextrecord();
	const std::string& s = A.currentrecordstring();
	//A.readnextgroup()
	std::vector<std::vector<int>> intfields;
	std::vector<std::vector<double>> doublefields;
	cStopWatch sw;
	size_t i1 = A.readnextgroup(4, intfields, doublefields);
	size_t i2 = A.readnextgroup(4, intfields, doublefields);
	size_t i3 = A.readnextgroup(4, intfields, doublefields);
	size_t i4 = A.readnextgroup(4, intfields, doublefields);
	sw.reportnow();
	prompttocontinue();
	return true;
};

bool test_aseggdfheader(){			
	std::string dfnpath = R"(z:\projects\earth_sci_test\test_data\output\inversion.output.dfn)";	
	cASEGGDF2Header H(dfnpath);
	H.write(dfnpath + ".txt");
	return true;
};

void test_marray(){	
	std::vector<size_t> dims = { 3, 4, 2 };
	Marray<int> a(dims.data(), dims.data() + dims.size());
	
	for (size_t j = 0; j<a.size(); ++j) a(j) = j;
	
	
	std::cout << a.asString() << std::endl;
	std::cout << a.size() << std::endl;

	dims = { 3, 2, 4 };
	a.reshape(dims.data(), dims.data() + dims.size());

	std::cout << a.asString() << std::endl;
	std::cout << a.size() << std::endl;
};

void test_subsample() {
	std::string inncpath = R"(D:\inversion_netcdf\line_data_em\902467_1_Field_Survey_Data_20200826.nc)";
	std::string outncpath = R"(D:\inversion_netcdf\line_data_em\test.nc)";
	cGeophysicsNcFile infile(inncpath, NcFile::FileMode::read);
	cGeophysicsNcFile outfile(outncpath, NcFile::FileMode::replace);

	std::vector<std::string> include_varnames = { "proj_client", "flight", "flight_index", "longitude", "latitude", "emz_nonhprg" };
	std::vector<std::string> exclude_varnames = { "easting", "northing" };
	outfile.subsample(infile,30,include_varnames, exclude_varnames);
}

void test_mpi(int argc, char** argv) {
	
	int mpisize, mpirank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
		
	//std::string inncpath = R"(D:\inversion_netcdf\line_data_em\902467_1_Field_Survey_Data_20200826.nc)";	
	std::string inncpath = R"(D:\inversion_netcdf\line_data_em\test.nc)";	
	
	if (mpirank == 0) {
		cGeophysicsNcFile nc(inncpath, NcFile::FileMode::write);
		nc.addSampleVar("newvar", NcType::nc_FLOAT);
		nc.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	cGeophysicsNcFile nc(inncpath, NcFile::FileMode::write);		
	cSampleVar var = nc.getSampleVar("newvar");
	for (size_t li = 0; li < nc.nlines(); li++) {
		if ((li%mpisize) == mpirank){
			float v = (mpirank + 1) * 10000 + li;
			std::vector<float> vals(nc.nlinesamples(li), v);
			var.putLine(li, vals);
			std::cout << "Rank " << mpirank << " Lineindex " << li << " v=" << v << std::endl;						
		}		
	}
	nc.sync();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

int main(int argc, char** argv)
{
	_GSTITEM_
	glog.logmsg("Opening log file\n");
	glog.open("test.log");	
	try{	
		test_mpi(argc,argv);
		//test_subsample();
		//example_magnetics();
		//example_aem_conductivity();	
		//test_create();
		//test_update();
		//test_aseggdfexport_1d();
		//test_aseggdfexport_2d();
		//test_columnfile();
		//test_aseggdfheader();
		//test_marray();
		//test_convert();		
		glog.close();
	}
	catch (NcException& e)
	{
		_GSTPRINT_		
		glog.logmsg(e.what());
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		glog.logmsg(e.what());
		return 1;
	}

	return 0;
}


