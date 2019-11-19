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

using namespace netCDF;
using namespace netCDF::exceptions;

#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
class cStackTrace gtrace;
#endif

#define _PROGRAM_ "intrepid2netcdf"
#define _VERSION_ "1.0"

#include <mpi.h>

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "intrepid.h"

#include "metadata.h"
#include "crs.h"
#include "csvfile.h"
#include "logger.h"
#include "geophysics_netcdf.h"

class cLogger glog; //The instance of the global log file manager

/*
class cIntrepidConverter {

	int mpisize;
	int mpirank;
	std::string IDBPath;
	std::string IDBName;
	std::string NCPath;

public:
	std::string LogFile;
	cMetaDataTable  A;
	cMetaDataRecord M;
	cCSVFile V;
	std::string IntrepidDatbasesDir;
	std::string NetCDFOutputDir;
	bool OverWriteExistingNcFiles = false;
	bool DummyRun = false;
	bool AddGeospatialMetadata = true;
	bool AddLineStartEndPoints = true;
	bool AddAlphaShapePolygon = true;

	size_t DatabaseStart = 0;
	size_t DatabaseSubsample = 1;

	cIntrepidConverter(const std::string& controlfile, int _mpisize = 1, int _mpirank = 0) {
		_GSTITEM_

		mpisize = _mpisize;
		mpirank = _mpirank;
		cBlock B(controlfile);
		LogFile = B.getstringvalue("LogFile");

		std::string suffix = stringvalue(mpirank, ".%04lu");
		LogFile = insert_after_filename(LogFile, suffix);

		glog.open(LogFile);
		glog.logmsg("Log file opened at %s\n", timestamp().c_str());
		glog.logmsg("Program %s \n", _PROGRAM_);
		glog.logmsg("Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg("Control file %s\n", controlfile.c_str());
		glog.logmsg("Processes %zu\n", mpisize);
		glog.logmsg("Rank      %zu\n", mpirank);
		glog.log(B.get_as_string());

		M = cMetaDataRecord(B, "Metadata");
		std::string argusfile = B.getstringvalue("ArgusMetaData");
		A = cMetaDataTable(argusfile, 1);
		V = cCSVFile(B.getstringvalue("VocabularyFile"));
		DatabaseStart = B.getsizetvalue("DatabaseStart");
		DatabaseSubsample = B.getsizetvalue("DatabaseSubsample");
		OverWriteExistingNcFiles = B.getboolvalue("OverWriteExistingNcFiles");
		DummyRun = B.getboolvalue("DummyRun");
		AddGeospatialMetadata = B.getboolvalue("AddGeospatialMetadata");
		AddLineStartEndPoints = B.getboolvalue("AddLineStartEndPoints");
		AddAlphaShapePolygon = B.getboolvalue("AddAlphaShapePolygon");
		IntrepidDatbasesDir = B.getstringvalue("IntrepidDatabasesPath");
		NetCDFOutputDir = B.getstringvalue("NetCDFOutputDir");
		fixseparator(IntrepidDatbasesDir);
		fixseparator(NetCDFOutputDir);

		process_databases();

		glog.logmsg("Finished at %s\n", timestamp().c_str());
		glog.close();
	};

	~cIntrepidConverter() {


	};

	bool process_databases() {

		std::vector<std::string> dlist = getfilelist(IntrepidDatbasesDir, ".DIR");
		int count = -1;
		for (size_t di = DatabaseStart - 1; di < dlist.size(); di = di + DatabaseSubsample) {
			count++;
			if (count%mpisize == mpirank) {
				IDBPath = ILDataset::dbdirpath(dlist[di]);
				IDBName = ILDataset::dbname(dlist[di]);
				NCPath = NetCDFOutputDir + IDBName + ".nc";

				if (exists(NCPath)) {
					if (OverWriteExistingNcFiles == false) {
						glog.logmsg("Warning 0: NetCDF file %s already exists - skipping this database\n", NCPath.c_str());
						continue;
					}
				}

				glog.logmsg("\nConverting database %zu of %zu: %s\n", di + 1, dlist.size(), dlist[di].c_str());

				bool datasetexists = exists(IDBPath);
				if (datasetexists == false) {
					glog.logmsg("Error 0: Database %zu does not exist: %s\n", di, IDBPath.c_str());
					continue;
				}

				glog.logmsg("Opening Intrepid database\n");
				ILDataset D(IDBPath);
				if (D.valid == false) {
					glog.logmsg("Error 1: problem opening database = skipping this database\n");
					continue;
				}

				if (D.fieldexists("longitude_GDA94") == false) {
					glog.logmsg("Error 2: field longitude_GDA94 does not exist = skipping this database\n");
					continue;
				}

				if (D.fieldexists("latitude_GDA94") == false) {
					glog.logmsg("Error 3: field latitude_GDA94 does not exist = skipping this database\n");
					continue;
				}

				if (D.hassurveyinfoid_fieldexists("X") == false) {
					glog.logmsg("Warning 1: could not determine the X field in the SurveyInfo file\n");
				}

				if (D.hassurveyinfoid_fieldexists("Y") == false) {
					glog.logmsg("Warning 2: could not determine the Y field in the SurveyInfo file\n");
				}

				glog.logmsg("Determining GA project number\n");
				int projectnumber = determineprojectnumber(D);
				if (projectnumber < 0) {
					glog.logmsg("Error 3: could not determinine GA project number - skipping this database\n");
					continue;
				}

				glog.logmsg("Populating metadata\n");
				cMetaDataRecord R = M;
				bool status = populate_metadata(R, projectnumber);
				if (status == false) {
					glog.logmsg("Error 4: could not find the Argus metadata record for this project %d - skipping this database\n", projectnumber);
					continue;
				}

				glog.logmsg("Creating NetCDF file\n");
				cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

				glog.logmsg("Adding line index variables\n");
				bool lstatus = add_lineindex(ncFile, D);
				if (lstatus == false) {
					glog.logmsg("Error 5: could not determine the line numbers\n");
					continue;
				}

				glog.logmsg("Adding global metadata\n");
				add_global_metadata(ncFile, R);

				glog.logmsg("Adding groupbyline varaibles\n");
				add_groupbyline_variables(ncFile, D);

				glog.logmsg("Adding indexed varaibles\n");
				add_indexed_variables(ncFile, D);

				if (DummyRun == false) {
					if (AddGeospatialMetadata) {
						glog.logmsg("Adding geospatial metadata\n");
						add_geospatial_metadata(ncFile, D);
					}

					if (AddLineStartEndPoints) {
						glog.logmsg("Adding line start and end points\n");
						ncFile.addLineStartEndPointsLL();
					}

					if (AddAlphaShapePolygon) {
						glog.logmsg("Adding alphashape polygon\n");
						ncFile.addAlphaShapePolygon("longitude", "latitude");
					}
				}
				glog.logmsg("Conversion complete\n");
			}
		}
		return true;
	}

	int determineprojectnumber(ILDataset& D) {

		_GSTITEM_
			bool hassurvey = false;
		std::string sfname = D.fieldnamelike("survey");
		if (sfname.size() == 0) {
			sfname = D.fieldnamelike("project");
		}

		if (sfname.size() == 0) {
			sfname = D.fieldnamelike("geoscience_australia_project_number");
		}

		if (sfname.size() > 0) {
			hassurvey = true;
		}

		int pmin = 0;
		if (hassurvey) {
			cStats<double> s = D.fieldstats(sfname);
			pmin = (int)s.min;
			if (s.min != s.max) {
				glog.logmsg("Warning 1a: %s Database has more than one survey number (%d to %d)\n", IDBPath.c_str(), (int)s.min, (int)s.max);
			}
		}

		//Check for GA or GSV project number in database name
		int  GApnumberguess = -1;
		bool GAguessed = parseGAprojectnumber(IDBName, GApnumberguess);
		if (GAguessed == false) {
			int  GSVpnumberguess = -1;
			bool GSVguessed = parseGSVprojectnumber(IDBName, GSVpnumberguess);
			if (GSVguessed) {
				size_t sindex = (size_t)A.findkeyindex("SURVEYNAME");
				for (size_t i = 0; i < A.records.size(); i++) {
					std::string sname = A.records[i][sindex];
					int gsvnum;
					bool status = parseGSVprojectnumber(sname, gsvnum);
					if (status) {
						if (gsvnum == GSVpnumberguess) {
							size_t pindex = (size_t)A.findkeyindex("PROJECT");
							GApnumberguess = std::atoi((A.records[i][pindex]).c_str());
							GAguessed = true;
							glog.logmsg("Warning 1b: GA Project number taken from GSV number (GSV=%d -> GA=%d)\n", GSVpnumberguess, GApnumberguess);
							break;
						}
					}
				}
			}
		}

		if (GAguessed == true) {
			if (hassurvey == true && pmin != GApnumberguess) {
				glog.logmsg("Warning 1c: Project number guessed from database name does not match that from the survey field in database (%d and %d) using %d\n", GApnumberguess, pmin, GApnumberguess);
			}
			return GApnumberguess;
		}
		else {
			if (hassurvey == false) {
				glog.logmsg("Warning 1d: Could not determine a project from the database name or survey field\n");
				return -1;
			}
			return pmin;
		}
	}

	bool parseGAprojectnumber(const std::string& str, int& pguess) {

		pguess = -1;
		for (size_t i = 0; i < IDBName.size(); i++) {
			if (sscanf(&str[i], "P%d", &pguess) == 1) {
				return true;
			}
			else if (sscanf(&str[i], "p%d", &pguess) == 1) {
				return true;
			}
		}
		return false;

	}

	bool parseGSVprojectnumber(const std::string& str, int& pguess) {
		pguess = -1;
		for (size_t i = 0; i < str.size(); i++) {
			if (sscanf(&str[i], "GSV%d", &pguess) == 1) {
				return true;
			}
			else if (sscanf(&str[i], "gsv%d", &pguess) == 1) {
				return true;
			}
		}
		return false;
	}

	bool populate_metadata(cMetaDataRecord& r, const int projectnumber) {
		_GSTITEM_
			std::vector<size_t> arecords = A.findmatchingrecords("PROJECT", projectnumber);

		if (arecords.size() == 0) {
			glog.logmsg("Warning 2a: Could not find a matching PROJECT for project number %d\n", projectnumber);
			return false;
		}

		if (arecords.size() > 1) {
			glog.logmsg("Warning 2b: Found more than 1 matching PROJECT for project number %d\n", projectnumber);
			return false;
		}
		size_t argusrecord = arecords[0];

		for (size_t i = 0; i < r.header.size(); i++) {
			std::string& s = r.values[i];
			if (std::strncmp(s.c_str(), "Argus.", 6) == 0) {
				std::string argusfield = s.substr(6, s.size() - 6);
				int ki = A.findkeyindex(argusfield);
				s = A.records[argusrecord][(size_t)ki];
			}
		}

		int ki = r.findkeyindex("date_created");
		if (ki >= 0) {
			std::string datestr = timestring("%Y-%m-%d %H:%M:%S");
			r.values[(size_t)ki] = datestr;
		}

		ki = r.findkeyindex("geoscience_australia_source_dataset");
		if (ki >= 0) {
			r.values[(size_t)ki] = IDBName;
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
			else {
				std::string msg = "Error 6: Unknown Intrepid data type\n";
				glog.logmsg(msg.c_str());
				throw(msg);
			}
	}

	bool set_intrepid_nullvalue(cGeophysicsVar& v)
	{
		_GSTITEM_
			NcType t = v.getType();
		if (t == ncUbyte) v.add_missing_value(IDatatype::bytenull());
		else if (t == ncShort) v.add_missing_value(IDatatype::shortnull());
		else if (t == ncInt) v.add_missing_value(IDatatype::intnull());
		else if (t == ncFloat) v.add_missing_value(v.preferred_float_missing_value());
		else if (t == ncDouble) v.add_missing_value(v.preferred_double_missing_value());
		else {
			std::string msg = "Error 7: Unsupported data type\n";
			glog.logmsg(msg.c_str());
			throw(msg);
		}
		return true;
	}

	bool add_lineindex(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
			std::vector<size_t> count = D.linesamplecount();
		std::vector<size_t> linenumbers;
		bool status = D.getlinenumbers(linenumbers);
		if (status) {
			ncFile.InitialiseNew(linenumbers, count);
			return true;
		}
		return false;
	}

	bool add_groupbyline_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
			if (D.valid == false)return false;
		size_t nlines = D.nlines();
		NcDim  dim_line = ncFile.getDim("line");

		size_t i_field = (size_t)V.findkeyindex("field_name");
		size_t i_convert = (size_t)V.findkeyindex("convert");
		size_t i_variable_name = (size_t)V.findkeyindex("variable_name");
		size_t i_standard_name = (size_t)V.findkeyindex("standard_name");
		size_t i_units = (size_t)V.findkeyindex("units");

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == false) continue;

			if (F.datatype().name() == "UNKNOWN") {
				glog.logmsg("Warning 5: skipping field %s: unsupported Intrepid datatype\n", F.Name.c_str());
				continue;
			}

			std::vector<size_t> recs = V.findmatchingrecords(i_field, F.Name);
			if (recs.size() == 0) {
				glog.logmsg("Warning 6: skipping field %s: could not find a match in the vocabulary file\n", F.Name.c_str());
				continue;
			}

			size_t vi = recs[0];
			std::string convert = V.records[vi][i_convert];
			if (convert == "no") {
				glog.logmsg("Warning 7: skipping ignored field %s\n", F.Name.c_str());
				continue;
			}
			glog.logmsg("Converting field %s\n", F.Name.c_str());

			std::string variable_name = V.records[vi][i_variable_name];
			std::string standard_name = V.records[vi][i_standard_name];
			std::string units = V.records[vi][i_units];

			if (variable_name == DN_LINE)continue;

			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cLineVar var = ncFile.addLineVar(variable_name, nc_datatype(F), dims);
			set_intrepid_nullvalue(var);
			var.add_long_name(standard_name);
			var.add_original_name(F.Name);
			var.add_units(units);

			for (size_t li = 0; li < nlines; li++) {
				ILSegment& S = F.Segments[li];

				if (S.readbuffer() == false) {
					glog.logmsg("Error 8: could not read buffer for line sequence number %zu in field %s\n", li, F.Name.c_str());
					return false;
				}

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

		size_t i_field = (size_t)V.findkeyindex("field_name");
		size_t i_convert = (size_t)V.findkeyindex("convert");
		size_t i_variable_name = (size_t)V.findkeyindex("variable_name");
		size_t i_standard_name = (size_t)V.findkeyindex("standard_name");
		size_t i_units = (size_t)V.findkeyindex("units");

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == true)continue;

			if (F.datatype().name() == "UNKNOWN") {
				glog.logmsg("Warning 5: skipping field %s: unsupported Intrepid datatype\n", F.Name.c_str());
				continue;
			}

			std::vector<size_t> recs = V.findmatchingrecords(i_field, F.Name);
			if (recs.size() == 0) {
				glog.logmsg("Warning 6: skipping field %s: could not find a match in the vocabulary file\n", F.Name.c_str());
				continue;
			}

			size_t vi = recs[0];
			std::string convert = V.records[vi][i_convert];
			if (convert == "no") {
				glog.logmsg("Warning 7: skipping ignored field %s\n", F.Name.c_str());
				continue;
			}
			glog.logmsg("Converting field %s\n", F.Name.c_str());

			std::string variable_name = V.records[vi][i_variable_name];
			std::string standard_name = V.records[vi][i_standard_name];
			std::string units = V.records[vi][i_units];

			NcVar tmp = ncFile.getVar(variable_name);
			if (tmp.isNull() == false) {
				glog.logmsg("Error 9: variable name %s for field %s already exists in this NC file - skipping this field\n", variable_name.c_str(), F.Name.c_str());
				continue;
			}

			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cSampleVar var = ncFile.addSampleVar(variable_name, nc_datatype(F), dims);
			set_intrepid_nullvalue(var);
			var.add_long_name(standard_name);
			var.add_original_name(F.Name);
			var.add_units(units);

			if (DummyRun) continue;

			size_t startindex = 0;
			for (size_t li = 0; li < nlines; li++) {
				ILSegment& S = F.Segments[li];

				if (S.readbuffer() == false) {
					glog.logmsg("Error 8: could not read buffer for line sequence number %zu in field %s\n", li, F.Name.c_str());
					return false;
				}

				//Replace the nulls with more "sensible" values
				if (S.datatype().isfloat()) {
					S.change_nullvalue(var.preferred_float_missing_value());
				}
				else if (S.datatype().isdouble()) {
					S.change_nullvalue(var.preferred_double_missing_value());
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

	bool add_global_metadata(cGeophysicsNcFile& ncFile, const cMetaDataRecord& m) {
		_GSTITEM_
			for (size_t j = 0; j < m.header.size(); j++) {
				if (m.header[j][0] == '/')continue;

				if (strcasecmp(m.header[j], "geoscience_australia_airborne_survey_project_number") == 0) {
					int pnum = atoi(m.values[j].c_str());
					ncFile.putAtt(m.header[j].c_str(), ncInt, pnum);
				}
				else if (strcasecmp(m.header[j], "nominal_minimum_line_spacing") == 0) {
					std::string s = m.values[j] + " m";
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else if (strcasecmp(m.header[j], "nominal_maximum_line_spacing") == 0) {
					std::string s = m.values[j] + " m";
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else if (strcasecmp(m.header[j], "nominal_height_above_ground") == 0) {
					std::string s = m.values[j] + " m";
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else if (strcasecmp(m.header[j], "nominal_height_above_sealevel") == 0) {
					std::string s = m.values[j] + " m";
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else if (strcasecmp(m.header[j], "acquisition_start_date") == 0) {
					std::string s = standardize_date(m.values[j]);
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else if (strcasecmp(m.header[j], "acquisition_end_date") == 0) {
					std::string s = standardize_date(m.values[j]);
					ncFile.putAtt(m.header[j].c_str(), s);
				}
				else {
					ncFile.putAtt(m.header[j].c_str(), m.values[j].c_str());
				}
			}
		return true;
	}

	std::string standardize_date(const std::string& indate) {
		_GSTITEM_
			int day, month, year;
		int n = sscanf(indate.c_str(), "%d/%d/%d", &day, &month, &year);
		if (n == 3) {
			std::string s = strprint("%04d-%02d-%02d", year, month, day);
			return s;
		}
		else return indate;
	}

	bool add_geospatial_metadata(cGeophysicsNcFile& ncFile, ILDataset& D) {
		_GSTITEM_

			std::string datum;
		std::vector<std::string> n;
		n.resize(0);
		n.push_back("longitude_GDA94");
		n.push_back("longitude_WGS84");
		n.push_back("longitude_AGD84");
		n.push_back("longitude_AGD66");
		for (size_t i = 0; i < n.size(); i++) {
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size() > 0) {
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
		for (size_t i = 0; i < n.size(); i++) {
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size() > 0) {
				cStats<double> st = D.fieldstats(fn);
				ncFile.putAtt("geospatial_lat_min", ncDouble, st.min);
				ncFile.putAtt("geospatial_lat_max", ncDouble, st.max);
				ncFile.putAtt("geospatial_lat_units", "degree_north");
				ncFile.putAtt("geospatial_lat_resolution", "point");
				break;
			}
		}

		if (datum.size() > 0) {
			int epsgcode = getepsgcode(datum);
			ncFile.addCRS(epsgcode);
		}

		n.resize(0);
		n.push_back("altitude");
		n.push_back("clearance");
		for (size_t i = 0; i < n.size(); i++) {
			std::string fn = D.fieldnamelike(n[i]);
			if (fn.size() > 0) {
				cStats<double> st = D.fieldstats(fn);
				ncFile.putAtt("geospatial_vertical_min", ncDouble, st.min);
				ncFile.putAtt("geospatial_vertical_max", ncDouble, st.max);
				ncFile.putAtt("geospatial_vertical_units", "m");
				ncFile.putAtt("geospatial_vertical_resolution", "point");
				ncFile.putAtt("geospatial_vertical_positive", "up");
				break;
			}
		}
		return true;
	}
};
*/

class cIntrepidToNetCDFConverter {
	int MPISize;
	int MPIRank;
	std::string IntrepiDatabasePath;
	std::string NCPath;
	bool OverWriteExistingNcFiles = true;
	std::string LogFile;		

public:

	cIntrepidToNetCDFConverter(const std::string& intrepiddatabasepath, const std::string& ncfilepath, const int mpisize, const int mpirank) {
		_GSTITEM_
		MPISize = mpisize;
		MPIRank = mpirank;
		IntrepiDatabasePath = fixseparator(intrepiddatabasepath);
		NCPath = fixseparator(ncfilepath);
		LogFile = NCPath + ".log";

		glog.open(LogFile);		
		glog.logmsg("Program %s \n", _PROGRAM_);
		glog.logmsg("Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		bool status = process();	
		if (status == false) {
			glog.logmsg(MPIRank,"Error 0: converting %s to %s\n",intrepiddatabasepath.c_str(), ncfilepath.c_str());
		}
		glog.close();
	};

	~cIntrepidToNetCDFConverter() {


	};

	bool process() {

		std::string IDBPath = ILDataset::dbdirpath(IntrepiDatabasePath);
		std::string IDBName = ILDataset::dbname(IntrepiDatabasePath);

		if (exists(NCPath)){
			if (OverWriteExistingNcFiles == false) {
				glog.logmsg("Warning 1: NetCDF file %s already exists - skipping this database\n", NCPath.c_str());
				return false;
			}
		}

		glog.logmsg("\nConverting database: %s\n", IntrepiDatabasePath.c_str());
		bool datasetexists = exists(IntrepiDatabasePath);
		if (datasetexists == false) {
			glog.logmsg("Error 1: Database does not exist: %s\n", IntrepiDatabasePath.c_str());
			return false;
		}

		glog.logmsg("\nOpening Intrepid database\n");
		ILDataset D(IntrepiDatabasePath);
		if (D.ispointdataset()) {
			glog.logmsg("Error 1: point databases not supported - skipping %s\n",IntrepiDatabasePath.c_str());
			return false;
		}
		else if(D.valid == false){
			glog.logmsg("Error 2: problem opening database - skipping %s\n",IntrepiDatabasePath.c_str());
			return false;
		}

		std::string linenumberfieldname;
		D.getlinenumberfieldname(linenumberfieldname);
		if (D.fieldexists(linenumberfieldname) == false) {
			glog.logmsg("Error 2: could not determine the LineNumber field in the SurveyInfo file\n");
		}

		if (D.hassurveyinfoid_fieldexists("X") == false) {
			glog.logmsg("Warning 1: could not determine the X field in the SurveyInfo file\n");
		}

		if (D.hassurveyinfoid_fieldexists("Y") == false) {
			glog.logmsg("Warning 2: could not determine the Y field in the SurveyInfo file\n");
		}

		glog.logmsg("Creating NetCDF file: %s\n",NCPath.c_str());
		cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

		glog.logmsg("\nAdding the line index variable\n");
		bool lstatus = add_lineindex(ncFile, D);
		if (lstatus == false) {
			glog.logmsg("Error 3: could not determine the line numbers\n");
			return false;
		}

		glog.logmsg("\nAdding groupby varaibles\n");
		add_groupbyline_variables(ncFile, D);

		glog.logmsg("\nAdding indexed varaibles\n");
		add_indexed_variables(ncFile, D);

		glog.logmsg("\nConversion complete\n");
		return true;
	}

	NcType nc_datatype(const ILField& F)
	{
		_GSTITEM_
		if (F.datatype().isbyte()) return NcType(ncUbyte);
		else if (F.datatype().isshort()) return NcType(ncShort);
		else if (F.datatype().isint()) return NcType(ncInt);
		else if (F.datatype().isfloat())return NcType(ncFloat);
		else if (F.datatype().isdouble())return NcType(ncDouble);
		else {
			std::string msg = "Error 6: Unknown Intrepid data type\n";
			glog.logmsg(msg.c_str());
			throw(msg);
		}
	}

	bool set_intrepid_nullvalue(cGeophysicsVar& v)
	{
		_GSTITEM_
		NcType t = v.getType();
		if (t == ncUbyte) v.add_missing_value(IDatatype::bytenull());
		else if (t == ncShort) v.add_missing_value(IDatatype::shortnull());
		else if (t == ncInt) v.add_missing_value(IDatatype::intnull());
		//UINT is not an Intrepid datatype //
		//else if (t == ncUint) v.add_missing_value(NC_FILL_UINT);
		else if (t == ncFloat) v.add_missing_value(IDatatype::floatnull());
		else if (t == ncDouble) v.add_missing_value(IDatatype::doublenull());
		else {
			std::string msg = strprint("Error 7: Unsupported data type %s for field %s\n",t.getName().c_str(), v.getName().c_str());
			glog.logmsg(msg.c_str());
			throw(msg);
		}
		return true;
	}

	bool add_lineindex(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		std::vector<size_t> count = D.linesamplecount();
		std::vector<size_t> linenumbers;
		bool status = D.getlinenumbers(linenumbers);
		if (status) {
			ncFile.InitialiseNew(linenumbers, count);
			return true;
		}
		return false;
	}

	bool add_groupbyline_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();
		NcDim  dim_line = ncFile.getDim("line");

		std::string linenumberfield;
		D.getlinenumberfieldname(linenumberfield);
		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == false) continue;

			if (F.datatype().name() == "UNKNOWN") {
				glog.logmsg("Warning 5: skipping field %s: unsupported Intrepid datatype\n", F.Name.c_str());
				continue;
			}

			glog.logmsg("Converting field %s\n", F.Name.c_str());
			
			if (F.Name == linenumberfield) continue;

			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cLineVar var = ncFile.addLineVar(F.Name, nc_datatype(F), dims);
			//set_intrepid_nullvalue(var);

			for (size_t li = 0; li < nlines; li++) {
				ILSegment& S = F.Segments[li];

				if (S.readbuffer() == false) {
					glog.logmsg("Error 8: could not read buffer for line sequence number %zu in field %s\n", li, F.Name.c_str());
					return false;
				}				
				change_fillvalues(S);

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

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == true)continue;

			if (F.datatype().name() == "UNKNOWN") {
				glog.logmsg("Warning 5: skipping field %s: unsupported Intrepid datatype\n", F.Name.c_str());
				continue;
			}

			glog.logmsg("Converting field %s\n", F.Name.c_str());
			NcVar tmp = ncFile.getVar(F.Name);
			if (tmp.isNull() == false) {
				glog.logmsg("Error 9: variable name %s already exists in this NC file - skipping this field\n", F.Name.c_str());
				return false;
			}

			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			cSampleVar var = ncFile.addSampleVar(F.Name, nc_datatype(F), dims);
			//set_intrepid_nullvalue(var);

			size_t startindex = 0;
			for (size_t li = 0; li < nlines; li++) {
				ILSegment& S = F.Segments[li];

				if (S.readbuffer() == false) {
					glog.logmsg("Error 8: could not read buffer for line sequence number %zu in field %s\n", li, F.Name.c_str());
					return false;
				}
				change_fillvalues(S);
				
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
	
	void change_fillvalues(ILSegment& S) {
		//Replace the nulls with more "sensible" values
		if (S.datatype().isbyte()) {
			S.change_nullvalue(defaultmissingvalue(ncByte));
		}
		else if (S.datatype().isshort()) {
			S.change_nullvalue(defaultmissingvalue(ncShort));
		}
		else if (S.datatype().isint()) {
			S.change_nullvalue(defaultmissingvalue(ncInt));
		}
		else if (S.datatype().isfloat()) {
			S.change_nullvalue(defaultmissingvalue(ncFloat));
		}
		else if (S.datatype().isdouble()) {
			S.change_nullvalue(defaultmissingvalue(ncDouble));
		}
	}

	bool add_global_metadata() {
		_GSTITEM_
		//ncFile.putAtt(m.header[j].c_str(), ncInt, pnum);		
		return true;
	}	
};

int main(int argc, char** argv)
{
	_GSTITEM_
	
	int mpisize;
	int mpirank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

	try
	{
		std::string dbdir    = argv[1];
		std::string ncdir    = argv[2];
		std::string listfile = argv[3];
		std::ifstream file(listfile);
		addtrailingseparator(dbdir);
		addtrailingseparator(ncdir);
		int k = 0;
		while (file.eof() == false) {
			std::string db;	
			file >> db;
			db = trim(db);
			if (db.size() > 0 && db[0] != '#') {
				sFilePathParts fpp = getfilepathparts(db);
				std::string dbname = dbdir + fpp.directory + fpp.prefix;
				std::string ncname = ncdir + fpp.directory + fpp.prefix + ".nc";				
				if (k % mpisize == mpirank) {
					std::cout << "[" << mpirank << "] " << dbname << " " << ncname << std::endl;
					cIntrepidToNetCDFConverter C(dbname, ncname, mpisize, mpirank);
				}
				k++;
			}
		}		
	}
	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		MPI_Finalize();
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		MPI_Finalize();
		return 1;
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


