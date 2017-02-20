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
#include <algorithm>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#include "stacktrace.h"
#ifdef USEGLOBALSTACKTRACE	
	cStackTrace globalstacktrace;
#endif

#include <fstream>
#include <mpi.h>
#include "general_utils.h"
#include "file_utils.h"
#include "asciicolumnfile.h"
#include "gdal_utils.h"

#include "metadata.h"
#include "crs.h"
#include "csvfile.h"
#include "geophysics_netcdf.h"

FILE* global_log_file = NULL;

class cASEGGDF2Converter{

	int mpisize;
	int mpirank;
	std::string DatPath;
	std::string DatName;
	std::string DfnPath;	
	std::string NCPath;

public:
	std::string LogFile;	
	
	cMetaDataTable   Argus;
	cMetaDataRecord  M;
	cCSVFile V;
	std::string DataFileDir;
	std::string NetCDFOutputDir;
	bool OverWriteExistingNcFiles = false;
	bool DummyRun = false;
	bool AddGeospatialMetadata = true;
	bool AddLineStartEndPoints = true;
	bool AddAlphaShapePolygon = true;

	size_t DataFileStart = 0;
	size_t DataFileSubsample = 1;

	cASEGGDF2Converter(const std::string& controlfile, int _mpisize = 1, int _mpirank = 0){		
		mpisize = _mpisize;
		mpirank = _mpirank;
		cBlock B(controlfile);
		LogFile = B.getstringvalue("LogFile");

		std::string suffix = stringvalue(mpirank, ".%04lu");
		LogFile = insert_after_filename(LogFile, suffix);

		global_log_file = fileopen(LogFile, "w");
		
		logmsg("Program starting at %s\n", timestamp().c_str());
		B.write(global_log_file);

		M = cMetaDataRecord(B, "Metadata");
		std::string argusfile = B.getstringvalue("ArgusMetaData");
		Argus = cMetaDataTable(argusfile, 1);
		V = cCSVFile(B.getstringvalue("VocabularyFile"));
		DataFileStart = B.getsizetvalue("DataFileStart");
		DataFileSubsample = B.getsizetvalue("DataFileSubsample");
		OverWriteExistingNcFiles = B.getboolvalue("OverWriteExistingNcFiles");
		DummyRun = B.getboolvalue("DummyRun");
		AddGeospatialMetadata = B.getboolvalue("AddGeospatialMetadata");
		AddLineStartEndPoints = B.getboolvalue("AddLineStartEndPoints");
		AddAlphaShapePolygon = B.getboolvalue("AddAlphaShapePolygon");

		DataFileDir = B.getstringvalue("ASEGGDF2Dir");
		NetCDFOutputDir = B.getstringvalue("NetCDFOutputDir");
		fixseparator(DataFileDir);
		fixseparator(NetCDFOutputDir);

		process();
	};

	~cASEGGDF2Converter(){						
		logmsg("Finished at %s\n", timestamp().c_str());
		fclose(global_log_file);
		global_log_file == NULL;
	};

	bool process(){
		std::vector<std::string> dlist = getfilelist(DataFileDir, "dat");
		int count = -1;
		for (size_t di = DataFileStart - 1; di < dlist.size(); di = di + DataFileSubsample){
			count++;
			if (count % mpisize == mpirank){
				DatPath = dlist[di];
				DatName = extractfilename_noextension(DatPath);
				DfnPath = extractfiledirectory(DatPath) + DatName + ".dfn";								
				NCPath = NetCDFOutputDir + DatName + ".nc";				
				convert_aseggdf2_file();
			}
		}
		return true;
	}

	size_t scan_for_line_index(const cAsciiColumnFile& D, std::vector<unsigned int>& line_index_start, std::vector<unsigned int>& line_index_count, std::vector<unsigned int>& line_number)
	{
		size_t fi = D.fieldindexbyname("line");
		size_t i1 = D.fields[fi].startchar;
		size_t i2 = D.fields[fi].endchar;

		size_t n = 0;
		std::ifstream in(DatPath);
		std::string s;

		int width = i2 - i1 + 1;
		unsigned int lastline = 0;
		size_t nlines = 0;
		while (std::getline(in, s)){
			std::string t = s.substr(i1, width);
			unsigned int lnum = atoi(t.data());
			if (lnum != lastline){
				line_number.push_back(lnum);
				line_index_start.push_back(n);
				line_index_count.push_back(1);
				lastline = lnum;
			}
			else{
				line_index_count.back()++;
			}
			n++;
		}
		return n;
	}

	std::vector<bool>  scan_for_groupby_fields(const cAsciiColumnFile& D, const std::vector<unsigned int>& line_index_count)
	{					
		std::vector<bool> groupby(D.fields.size(),true);
		std::ifstream infile(DatPath);
		std::string s, t;			
		int nl = std::min((int)2,(int)line_index_count.size());
		for (size_t li = 0; li < nl; li++){						
			std::getline(infile, s);
			for (size_t si = 1; si < line_index_count[li]; si++){				
				std::getline(infile, t);
				for (size_t fi = 0; fi < D.fields.size(); fi++){					
					if (groupby[fi] == true){						
						size_t i1 = D.fields[fi].startchar;
						size_t i2 = D.fields[fi].endchar;
						size_t width = i2 - i1 + 1;
						std::string a = s.substr(i1, width);
						std::string b = t.substr(i1, width);
						if (a != b) groupby[fi] = false;						
					}
				}
			}
		}			
		return groupby;		
	}

	int determine_GA_project_number(const cAsciiColumnFile& D)
	{		
		int project = -1;
		size_t fi = D.nullfieldindex();		
		if (D.fieldindexbyname("project", fi)){}
		else if (D.fieldindexbyname("ga_project", fi)){}
		else if (D.fieldindexbyname("survey", fi)){}
		else if (parse_GA_projectnumber(DatName, project)){
			return project;
		}
		else return -1;
					
		size_t n = 0;
		std::ifstream in(DatPath);
		std::string s;
		size_t i1 = D.fields[fi].startchar;
		size_t i2 = D.fields[fi].endchar;
		int width = i2 - i1 + 1;				
		std::getline(in, s);
		std::string t = s.substr(i1, width);
		project = atoi(t.data());		
		return project;
	}

	bool parse_GA_projectnumber(const std::string& str, int& pguess){
		bool status = false;
		pguess = -1;
		for (size_t i = 0; i < DatName.size(); i++){
			if (sscanf(&str[i], "P%d", &pguess) == 1){
				return true;
			}
			else if (sscanf(&str[i], "p%d", &pguess) == 1){
				return true;
			}
		}
		return false;
	}
	
	bool populate_metadata(cMetaDataRecord& r, const int projectnumber){		
		std::vector<size_t> arecords = Argus.findmatchingrecords("PROJECT", projectnumber);

		if (arecords.size() == 0){
			logmsg("Warning 2a: Could not find a matching PROJECT for project number %d\n", projectnumber);
			return false;
		}

		if (arecords.size() > 1){
			logmsg( "Warning 2b: Found more than 1 matching PROJECT for project number %d\n", projectnumber);
			return false;
		}
		size_t argusrecord = arecords[0];

		for (size_t i = 0; i < r.header.size(); i++){
			std::string& s = r.values[i];
			if (std::strncmp(s.c_str(), "Argus.", 6) == 0){
				std::string argusfield = s.substr(6, s.size() - 6);
				int ki = Argus.findkeyindex(argusfield);
				s = Argus.records[argusrecord][(size_t)ki];
			}
		}

		int ki = r.findkeyindex("date_created");
		if (ki >= 0){
			std::string datestr = timestring("%Y-%m-%d %H:%M:%S");
			r.values[(size_t)ki] = datestr;
		}

		ki = r.findkeyindex("geoscience_australia_source_dataset");
		if (ki >= 0){
			r.values[(size_t)ki] = DatName;
		}
		return true;
	};

	bool convert_aseggdf2_file(){

		if (exists(DatPath) == false){
			logmsg( "Error 1: Data file %s does not exist\n", DatPath.c_str());
			return false;
		}

		if (exists(DatPath) == false){
			logmsg( "Error 2: DFN file %s does not exist\n", DfnPath.c_str());
			return false;
		}

		if (exists(NCPath)){
			if (OverWriteExistingNcFiles == false){
				logmsg( "Warning 1: NetCDF file %s already exists - skipping this data file\n", NCPath.c_str());
				return false;
			}
		}

		logmsg( "\nConverting data file %s\n", DatPath.c_str());

		logmsg( "Opening data file\n");
		cAsciiColumnFile D(DatPath);
		
		logmsg( "Parsing ASEGGDF2 header\n");
		D.parse_aseggdf2_header(DfnPath);

		logmsg( "Determining GA project number\n");
		int projectnumber = determine_GA_project_number(D);
		if (projectnumber < 0){
			logmsg( "Error 3: could not determinine GA project number - skipping this database\n");
			return false;
		}

		logmsg( "Populating metadata\n");
		cMetaDataRecord R = M;
		bool status = populate_metadata(R, projectnumber);
		if (status == false){
			logmsg( "Error 4: could not find the Argus metadata record for this project %d - skipping this database\n", projectnumber);
			return false;
		}

		std::vector<unsigned int> line_number;
		std::vector<unsigned int> line_index_start;
		std::vector<unsigned int> line_index_count;
		logmsg( "Scanning for line index\n");
		size_t npoints = scan_for_line_index(D,line_index_start, line_index_count, line_number);
		logmsg( "Total number of points is %lu\n", npoints);

		logmsg( "Scanning for groupby fields\n");
		std::vector<bool> isgroupby = scan_for_groupby_fields(D,line_index_count);

		logmsg( "Creating NetCDF file\n");
		cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

		logmsg( "Adding line index variables\n");
		ncFile.InitialiseNew(line_number, line_index_count);

		logmsg( "Adding global metadata\n");
		add_global_metadata(ncFile, R);
		
		//Pre process the fields
		for (size_t fi = 0; fi < D.fields.size(); fi++){
			std::string fieldname = D.fields[fi].name;
			if (fieldname == DN_LINE){
				logmsg("Skip processing field %3lu - %s (already in index)\n", fi, fieldname.c_str());
				continue;
			}
			else{
				logmsg("Pre processing field  %3lu - %s\n", fi, fieldname.c_str());
			}
			
			
			
			size_t nbands = D.fields[fi].nbands;

			std::vector<NcDim> vardims;
			if (nbands > 1){
				std::string dimname = "layers";
				NcDim dimband = ncFile.getDim(dimname);
				if (dimband.isNull()){
					std::vector<unsigned int> layers = increment((unsigned int)nbands, (unsigned int)0, (unsigned int)1);
					dimband = ncFile.addDimVar("layers", layers);					
				}
				bool status = dimband.isNull();
				vardims.push_back(dimband);
			}

			std::string varname = fieldname;
			nc_type vartype;

			if (D.fields[fi].isinteger()){
				vartype = NC_INT;
			}
			else{
				if (strncasecmp(fieldname.c_str(), "lat", 3) == 0)vartype = NC_DOUBLE;
				else if (strncasecmp(fieldname.c_str(), "lon", 3) == 0)vartype = NC_DOUBLE;
				else if (strncasecmp(fieldname.c_str(), "east", 4) == 0)vartype = NC_DOUBLE;
				else if (strncasecmp(fieldname.c_str(), "north", 5) == 0)vartype = NC_DOUBLE;
				else vartype = NC_FLOAT;					
			}

			if (isgroupby[fi]){
				cLineVar var = ncFile.addLineVar(varname, vartype, vardims);
				var.add_standard_name(fieldname);
				var.add_original_name(fieldname);
				var.add_units(D.fields[fi].units);
			}
			else{
				cSampleVar var = ncFile.addSampleVar(varname, vartype, vardims);
				var.add_standard_name(fieldname);
				var.add_original_name(fieldname);
				var.add_units(D.fields[fi].units);
			}			
		}

		size_t fi_line = D.fieldindexbyname("line");
		size_t fi_x    = D.fieldindexbyname("easting");
		size_t fi_y    = D.fieldindexbyname("northing");
		
		std::vector<std::vector<int>>    intfields;
		std::vector<std::vector<double>> dblfields;
		size_t lineindex = 0;
		logmsg("Processing lines\n");
		while (size_t nsamples = D.readnextgroup(fi_line, intfields, dblfields)){			
			for (size_t fi = 0; fi < D.fields.size(); fi++){								
				std::string fieldname = D.fields[fi].name;
				if (fieldname == DN_LINE) continue;

				size_t nbands = D.fields[fi].nbands;

				NcVar var = ncFile.getSampleVar(fieldname);
				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);								
				for (size_t bi = 0; bi < nbands; bi++){
					size_t nactive;
					if (isgroupby[fi]){
						nactive = 1;
						startp[0] = lineindex;
						startp[1] = bi;
						countp[0] = 1;
						countp[1] = 1;
					}
					else{
						nactive = nsamples;
						startp[0] = line_index_start[lineindex];
						startp[1] = bi;
						countp[0] = line_index_count[lineindex];
						countp[1] = 1;
					}
															
					if (D.fields[fi].isinteger()){
						std::vector<int> data(nactive);
						for (size_t si = 0; si < nactive; si++){							
							data[si] = intfields[fi][si*nbands + bi];
							if (data[si] == D.fields[fi].nullvalue){
								data[si] = NcIntNull;
							}
						}
						var.putVar(startp, countp, data.data());
					}
					else{
						std::vector<double> data(nsamples);
						for (size_t si = 0; si < nactive; si++){
							data[si] = dblfields[fi][si*nbands + bi];
							if (data[si] == D.fields[fi].nullvalue){
								data[si] = NcDoubleNull;
							}
						}
						var.putVar(startp, countp, data.data());
					}
				}				
			}			
			lineindex++;
		}
		D.closefile();

		if (DummyRun == false){			
			if (AddGeospatialMetadata){
				logmsg("Adding geospatial metadata\n");
				ncFile.addGeospatialMetadataXY();
				ncFile.addGeospatialMetadataVertical();
			}

			if (AddLineStartEndPoints){
				logmsg( "Adding line start and end points\n");
				ncFile.addLineStartEndPointsEN();
			}

			if (AddAlphaShapePolygon){
				logmsg( "Adding alphashape polygon\n");
				//ncFile.addAlphaShapePolygon("longitude","latitude");
				ncFile.addAlphaShapePolygon("easting", "northing");
			}
		}
		logmsg( "Conversion complete\n");

		return true;
	}
	
	bool add_global_metadata(cGeophysicsNcFile& ncFile, const cMetaDataRecord& m){		
		for (size_t j = 0; j < m.header.size(); j++){
			if (m.header[j][0] == '/')continue;

			if (strcasecmp(m.header[j], "geoscience_australia_airborne_survey_project_number") == 0){
				int pnum = atoi(m.values[j].c_str());
				ncFile.putAtt(m.header[j].c_str(), ncInt, pnum);
			}
			else if (strcasecmp(m.header[j], "nominal_minimum_line_spacing") == 0){
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

	std::string standardize_date(const std::string& indate){		
		int day, month, year;
		int n = sscanf(indate.c_str(), "%d/%d/%d", &day, &month, &year);
		if (n == 3){
			std::string s = strprint("%04d-%02d-%02d", year, month, day);
			return s;
		}
		else return indate;
	}

};

int main(int argc, char** argv)
{
	_GSTITEM_
	if (argc != 2){
		printf("Usage: %s control_file_name\n", argv[0]);
		return 1;
	}

	int mpisize;
	int mpirank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

	try
	{
		double t1 = gettime();
		std::string controlfile = argv[1];
		cASEGGDF2Converter C(controlfile, mpisize, mpirank);		
		double t2 = gettime();
		printf("Reading data elapsed time = %.2lf\n", t2 - t1);		
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

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}





