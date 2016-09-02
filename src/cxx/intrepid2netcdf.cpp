/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <netcdf>
#include <vector>
#include <limits>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#ifdef USEGLOBALSTACKTRACE
	#include "stacktrace.h"
	cStackTrace globalstacktrace;
#endif

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "intrepid.h"

#include "metadata.h"
#include "crs.h"
#include "csvfile.h"

#include "geophysics_netcdf.h"

class cIntrepidConverter{
	
	std::string IDBPath;
	std::string IDBName;
	std::string NCPath;

public:
	std::string LogFile;
	FILE* flog;		
	cMetaDataTable  A;	
	cMetaDataRecord M;
	cCSVFile V;	
	std::string IntrepidDatbasesDir;	
	std::string NetCDFOutputDir;
	bool DummyRun;
	size_t DatabaseStart = 0;
	size_t DatabaseSubsample = 1;

	cIntrepidConverter(const std::string& controlfile){
		_GSTITEM_

		cBlock B(controlfile);
		LogFile = B.getstringvalue("LogFile");
		flog = fileopen(LogFile, "w");
		message(flog, "Program starting at %s\n", timestamp().c_str());
		B.write(flog);

		M = cMetaDataRecord(B, "Metadata");
		std::string argusfile = B.getstringvalue("ArgusMetaData");
		A = cMetaDataTable(argusfile, 1);
		V = cCSVFile(B.getstringvalue("VocabularyFile"));
		DummyRun = B.getboolvalue("DummyRun");
		DatabaseStart     = B.getsizetvalue("DatabaseStart");
		DatabaseSubsample = B.getsizetvalue("DatabaseSubsample");
		
		IntrepidDatbasesDir = B.getstringvalue("IntrepidDatabasesPath");
		NetCDFOutputDir     = B.getstringvalue("NetCDFOutputDir");
		fixseparator(IntrepidDatbasesDir);
		fixseparator(NetCDFOutputDir);
	};

	~cIntrepidConverter(){
		_GSTITEM_
		rootmessage("Finished\n");
		rootmessage(flog, "Finished at %s\n", timestamp().c_str());
		fclose(flog);
	};

	bool process_databases(){

		std::vector<std::string> dlist = getfilelist(IntrepidDatbasesDir, ".DIR");
		for (size_t di = DatabaseStart - 1; di < dlist.size(); di = di + DatabaseSubsample){
			IDBPath = ILDataset::dbdirpath(dlist[di]);
			IDBName = ILDataset::dbname(dlist[di]);
			NCPath  = NetCDFOutputDir + IDBName + ".nc";

			message(flog, "\nConverting database %lu of %lu: %s\n", di+1, dlist.size(), dlist[di].c_str());
			
			bool datasetexists = exists(IDBPath);
			if (datasetexists == false){
				message(flog, "Warning: Database %lu does not exist: %s\n", di, IDBPath.c_str());
				continue;
			}

			ILDataset D(IDBPath);

			int projectnumber = guessprojectnumber(D);
			if (projectnumber < 0) continue;

			cMetaDataRecord R = M;
			populate_metadata(R,projectnumber);
						
			cGeophysicsNcFile ncFile(NCPath, NcFile::replace);
			add_lineindex(ncFile, D);
			add_groupbyline_variables(ncFile, D);			
			add_indexed_variables(ncFile, D);
			add_line_startend_points(ncFile, D);			
			add_global_metadata(ncFile, R);
			//add_geospatial_metadata(ncFile, D);
			ncFile.addAlphaShapePolygon();
		}
		return true;
	}

	int guessprojectnumber(ILDataset& D){
		
		_GSTITEM_
		bool hassurvey = false;
		std::string sfname = D.fieldnamelike("survey");
		if (sfname.size() == 0){
			sfname = D.fieldnamelike("project");
		}

		if (sfname.size() == 0){
			sfname = D.fieldnamelike("geoscience_australia_project_number");
		}
		
		if (sfname.size() > 0){
			hassurvey = true;
		}
		
		int pmin=-1, pmax=-1;
		if (hassurvey){			
			cStats<double> s = D.fieldstats(sfname);
			pmin = (int)s.min;
			pmax = (int)s.max;
			if (s.min != s.max){
				message(flog, "Warning 1a: %s Database has more than one survey number (%d to %d)\n", IDBPath.c_str(), (int)s.min, (int)s.max);
			}			
		}

		int pguess = 0;
		bool guessed = false;
		for (size_t i = 0; i < IDBName.size(); i++){
			if (IDBName[i] == 'P'){
				if (sscanf(&IDBName[i], "P%d", &pguess) == 1){
					guessed = true;
					break;
				}
			}
			else if (IDBName[i] == 'p'){
				if (sscanf(&IDBName[i], "p%d", &pguess) == 1){
					guessed = true;
					break;
				}
			}
		}

		if(guessed == true){
			if(hassurvey==true && pmin!=pguess){
				message(flog, "Warning 1b: %s Project number guessed from database name does not match that from the survey field in database (%d and %d) using %d\n",IDBPath.c_str(),pguess,pmin,pguess);
			}
			return pguess;
		}
		else{
			if (hassurvey==false){
				message(flog, "Warning 1c: %s Could not determine a project from the database name or survey field\n", IDBPath.c_str());
				return -1;
			}			
			return pmin;
		}		
	}

	bool populate_metadata(cMetaDataRecord& r, const int projectnumber){
		_GSTITEM_
		std::vector<size_t> arecords = A.findmatchingrecords("PROJECT", projectnumber);

		if (arecords.size() == 0){
			message(flog, "Warning 2a: Could not find a matching PROJECT for project number %d\n", projectnumber);
			return false;
		}

		if (arecords.size() > 1){
			message(flog, "Warning 2b: Found more than 1 matching PROJECT for project number %d\n",projectnumber);
			return false;
		}
		size_t argusrecord = arecords[0];
		
		for (size_t i = 0; i < r.header.size(); i++){
			std::string& s = r.values[i];
			if (std::strncmp(s.c_str(), "Argus.", 6) == 0){
				std::string argusfield = s.substr(6, s.size() - 6);
				int ki = A.findkeyindex(argusfield);
				s = A.records[argusrecord][(size_t)ki];
			}
		}

		int ki = r.findkeyindex("date_created");
		if (ki >= 0){
			std::string datestr = timestring("%Y-%m-%d %H:%M:%S");
			r.values[(size_t)ki] = datestr;
		}

		ki = r.findkeyindex("geoscience_australia_source_dataset");
		if (ki >= 0){
			r.values[(size_t)ki] = IDBPath;
		}		
		return true;
	};

	NcType nc_datatype(const ILField& F)
	{
		_GSTITEM_
		if (F.datatype().isbyte()) return NcType(ncUbyte);
		else if (F.datatype().isshort()) return NcType(ncShort);
		else if (F.datatype().isint()) return NcType(ncInt);			
		else if (F.datatype().isfloat())return NcType(ncFloat);
		else if (F.datatype().isdouble())return NcType(ncDouble);			
		else{
			std::string msg = "Error: Unknown Intrepid data type\n";
			message(flog, msg.c_str());
			throw(msg);
		}
	}

	bool set_intrepid_nullvalue(cGeophysicsVar& v)
	{
		_GSTITEM_
		NcType t = v.getType();
		if (t == ncUbyte) v.add_fillvalue(IDatatype::bytenull());
		else if (t == ncShort) v.add_fillvalue(IDatatype::shortnull());
		else if (t == ncInt) v.add_fillvalue(IDatatype::intnull());		
		else if (t == ncFloat) v.add_fillvalue(v.preferred_float_fillvalue());
		else if (t == ncDouble) v.add_fillvalue(v.preferred_double_fillvalue());
		else{
			std::string msg = "Error: Unsupported data type\n";
			message(flog, msg.c_str());
			throw(msg);
		}
		return true;
	}	

	bool add_lineindex(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_						
		size_t nsamples = D.nsamples();
		size_t nlines   = D.nlines();
		std::vector<size_t> count = D.linesamplecount();
		std::vector<size_t> linenumbers;
		D.getlinenumbers(linenumbers);
		ncFile.InitialiseNew(linenumbers,count);		
		return true;
	}
	
	bool add_groupbyline_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();
		NcDim  dim_line = ncFile.getDim("line");

		int i_field = V.findkeyindex("field_name");
		int i_convert = V.findkeyindex("convert");
		int i_variable_name = V.findkeyindex("variable_name");
		int i_standard_name = V.findkeyindex("standard_name");
		int i_units = V.findkeyindex("units");

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it){
			fi++;
			ILField& F = *it;
			//message(flog,"Processing field %lu %s %s\n", fi, F.Name.c_str(), F.datatype().name().c_str());
			if (F.isgroupbyline() == false) continue;

			std::vector<size_t> recs = V.findmatchingrecords(i_field, F.Name);
			if (recs.size() == 0){
				message(flog, "Warning: skipping field %s: could not find a match in the vocabulary file\n", F.Name.c_str());
				continue;
			}

			size_t vi = recs[0];
			std::string convert = V.records[vi][i_convert];
			if (convert == "no"){
				message(flog, "Ignoring field %s\n", F.Name.c_str());
				continue;
			}
			message(flog, "Converting field %s\n", F.Name.c_str());

			std::string variable_name = V.records[vi][i_variable_name];
			std::string standard_name = V.records[vi][i_standard_name];
			std::string units = V.records[vi][i_units];

			if (variable_name == DN_LINE)continue;

			std::vector<NcDim> dims;			
			if (F.nbands() > 1){
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cLineVar var = ncFile.addLineVar(variable_name,nc_datatype(F),dims);			
			set_intrepid_nullvalue(var);
			var.add_standard_name(standard_name);			
			var.add_original_name(F.Name);
			var.add_units(units);						
			
			for (size_t li = 0; li < nlines; li++){
				ILSegment& S = F.Segments[li];
				S.readbuffer();
				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);

				//Point dimension
				startp[0] = li;
				countp[0] = 1;

				//Band dimension
				startp[1] = 0;
				countp[1] = S.nbands();
				var.putVar(startp, countp, S.pvoid_groupby());
			}
		}
		return true;
	}
	
	bool add_indexed_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();		

		int i_field = V.findkeyindex("field_name");
		int i_convert = V.findkeyindex("convert");
		int i_variable_name = V.findkeyindex("variable_name");
		int i_standard_name = V.findkeyindex("standard_name");
		int i_units = V.findkeyindex("units");

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it){
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == true)continue;

			if (F.datatype().name() == "UNKNOWN"){
				message(flog, "Warning: skipping field %s: unsupported Intrepid datatype\n", F.Name.c_str());
				continue;
			}

			std::vector<size_t> recs = V.findmatchingrecords(i_field, F.Name);
			if (recs.size() == 0){
				message(flog, "Warning: skipping field %s: could not find a match in the vocabulary file\n", F.Name.c_str());
				continue;
			}

			size_t vi = recs[0];
			std::string convert = V.records[vi][i_convert];
			if (convert == "no"){
				message(flog, "Ignoring field %s\n", F.Name.c_str());
				continue;
			}
			message(flog, "Converting field %s\n", F.Name.c_str());

			std::string variable_name = V.records[vi][i_variable_name];
			std::string standard_name = V.records[vi][i_standard_name];
			std::string units = V.records[vi][i_units];

			std::vector<NcDim> dims;			
			if (F.nbands() > 1){
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cSampleVar var = ncFile.addSampleVar(variable_name, nc_datatype(F), dims);
			set_intrepid_nullvalue(var);
			var.add_standard_name(standard_name);
			var.add_original_name(F.Name);
			var.add_units(units);
						
			if (DummyRun) continue;

			size_t startindex = 0;
			for (size_t li = 0; li < nlines; li++){
				ILSegment& S = F.Segments[li];
				S.readbuffer();

				//Replace the nulls with more "sensible" values
				if (S.datatype().isfloat()){
					S.change_nullvalue(var.preferred_float_fillvalue());
				}
				else if (S.datatype().isdouble()){
					S.change_nullvalue(var.preferred_double_fillvalue());
				}

				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);

				//Sample dimension
				startp[0] = startindex;
				countp[0] = S.nsamples();

				//Band dimension
				startp[1] = 0;
				countp[1] = S.nbands();

				var.putVar(startp, countp, S.pvoid());
				startindex += S.nsamples();
			}
		}
		return true;
	}

	bool add_global_metadata(cGeophysicsNcFile& ncFile, const cMetaDataRecord& m){
		_GSTITEM_
		for (size_t j = 0; j < m.header.size(); j++){
			if (strcasecmp(m.header[j],"nominal_minimum_line_spacing")==0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_maximum_line_spacing") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_height_above_ground") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_height_above_sealevel") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "acquisition_start_date") == 0){
				std::string s = standardize_date(m.values[j]);
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "acquisition_end_date") == 0){
				std::string s = standardize_date(m.values[j]);
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else{
				ncFile.putAtt(m.header[j].c_str(), m.values[j].c_str());
			}			
		}		
		return true;
	}

	bool add_line_startend_points(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		std::vector<double> x1;
		std::vector<double> x2;
		std::vector<double> y1;
		std::vector<double> y2;
		D.get_line_start_end_points(x1,x2,y1,y2);
		ncFile.addLineStartEndPoints(x1, x2, y1, y2);		
		return true;
	}

	std::string standardize_date(const std::string& indate){
		_GSTITEM_
		int day, month, year;
		int n = sscanf(indate.c_str(), "%d/%d/%d", &day, &month, &year);
		if (n == 3){
			std::string s = strprint("%04d-%02d-%02d", year, month, day);			
			return s;
		}
		else return indate;
	}

	bool add_geospatial_metadata(cGeophysicsNcFile& ncFile, ILDataset& D){
		_GSTITEM_
		
		std::string datum;
		std::vector<std::string> n;
		n.resize(0);
		n.push_back("longitude_GDA94");
		n.push_back("longitude_WGS84");
		n.push_back("longitude_AGD84");
		n.push_back("longitude_AGD66");
		for (size_t i = 0; i < n.size(); i++){
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size()>0){
				cStats<double> st = D.fieldstats(fn);
				ncFile.putAtt("geospatial_lon_min", ncDouble, st.min);
				ncFile.putAtt("geospatial_lon_max", ncDouble, st.max);
				ncFile.putAtt("geospatial_lon_units", "degree_east");
				ncFile.putAtt("geospatial_lon_resolution", "point");
				if (i == 0)datum = "GDA94";
				else if (i == 0)datum = "WGS84";
				else if (i == 0)datum = "AGD84";
				else if (i == 0)datum = "AGD66";
				else datum = "";
				break;
			}
		}

		n.resize(0);
		n.push_back("latitude_GDA94");
		n.push_back("latitude_WGS84");
		n.push_back("latitude_AGD84");
		n.push_back("latitude_AGD66");
		for (size_t i = 0; i < n.size(); i++){
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size()>0){
				cStats<double> st = D.fieldstats(fn);
				ncFile.putAtt("geospatial_lat_min", ncDouble, st.min);
				ncFile.putAtt("geospatial_lat_max", ncDouble, st.max);
				ncFile.putAtt("geospatial_lat_units", "degree_north");
				ncFile.putAtt("geospatial_lat_resolution", "point");
				break;
			}
		}

		n.resize(0);
		n.push_back("altitude");
		n.push_back("clearance");
		for (size_t i = 0; i < n.size(); i++){
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size()>0){
				cStats<double> st = D.fieldstats(fn);
				ncFile.putAtt("geospatial_vertical_min", ncDouble, st.min);
				ncFile.putAtt("geospatial_vertical_max", ncDouble, st.max);
				ncFile.putAtt("geospatial_vertical_units", "m");
				ncFile.putAtt("geospatial_vertical_resolution", "point");
				ncFile.putAtt("geospatial_vertical_positive", "up");
				break;
			}
		}

		if (datum.size() > 0){
			cCRS crs(datum);	
			ncFile.addCRS(crs);			
		}

		return true;
	}	
		
};

int main(int argc, char** argv)
{	
	_GSTITEM_
	if (argc != 2){
		printf("Usage: %s control_file_name\n", argv[0]);
		return 1;
	}

	try
	{
		std::string controlfile = argv[1];
		cIntrepidConverter C(controlfile);
		C.process_databases();			
	}

	catch (NcException& e)
	{
		_GSTPRINT_		
		std::cout << e.what() << endl;
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << endl;
		return 1;
	}
	std::printf("Success\n");
	return 0;  // successfully terminated
}


