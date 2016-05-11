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

class cCSVFile{

	public:
		std::vector<std::string> header;
		std::vector<std::vector<std::string>> records;

		cCSVFile(){ }

		cCSVFile(const std::string csvfile){

			FILE* fp = fileopen(csvfile, "r");

			std::string str;
			std::vector<std::string> tokens;

			//Read header line
			filegetline(fp, str);
			split(str, ',', tokens);
			header = tokens;
			size_t nfields = header.size();			

			//Read remaining lines
			size_t k = 1;
			while (filegetline(fp, str)){
				k++;
				tokens.resize(0);
				split(str, ',', tokens);
				if (tokens.size() == nfields - 1){
					tokens.push_back("");
				}

				if (tokens.size() != nfields){
					std::printf("Error: On line %lu of file %s\n", k, csvfile.c_str());
					std::printf("Error: The number of header items (%lu) does not match the number of data items (%lu)\n", nfields, tokens.size());
				}

				for (size_t i = 0; i < tokens.size(); i++){
					tokens[i] = trim(tokens[i]);
					tokens[i] = stripquotes(tokens[i]);
				}
				records.push_back(tokens);
			}
			fclose(fp);
		}

		bool addfield(const std::string fname){
			header.push_back(fname);
			for (size_t i = 0; i < records.size(); i++){
				const std::string empty;
				records[i].push_back(empty);
			}
			return true;
		}

		bool setfield(const std::string fname, const std::string value, const size_t recindex){
			int k = findkeyindex(fname);
			if (k < 0)return false;
			records[recindex][(size_t)k] = value;
			return true;
		}

		bool setfield(const std::string fname, const std::string value){
			int k = findkeyindex(fname);
			if (k < 0)return false;
			if (records.size() == 0){
				records.resize(1);
				records[0].resize(header.size());
			}

			for (size_t i = 0; i < records.size(); i++){
				records[i][(size_t)k] = value;
			}
			return true;
		}

		int findkeyindex(const std::string& fname){
			for (size_t i = 0; i < header.size(); i++){
				if (header[i] == fname){
					return (int)i;
				}
			}
			return -1;
		}

		std::vector<size_t> findmatchingrecords(const size_t keyindex, const int value){
			std::vector<size_t> indices;
			for (size_t i = 0; i < records.size(); i++){
				int n = atoi(records[i][keyindex].c_str());
				if (n == value){
					indices.push_back(i);
				}
			}
			return indices;
		}

		std::vector<size_t> findmatchingrecords(const std::string key, const int value)	{
			int keyindex = findkeyindex(key);
			if (keyindex < 0){
				std::vector<size_t> indices;
				return indices;
			}
			return findmatchingrecords((size_t)keyindex, value);
		}

		std::vector<size_t> findmatchingrecords(const size_t keyindex, const std::string value){
			std::vector<size_t> indices;
			for (size_t i = 0; i < records.size(); i++){
				if (!strcasecmp(records[i][keyindex],value)){
					indices.push_back(i);
				}
			}
			return indices;
		}

		std::vector<size_t> findmatchingrecords(const std::string key, const std::string value)	{
			int keyindex = findkeyindex(key);
			if (keyindex < 0){
				std::vector<size_t> indices;
				return indices;
			}
			return findmatchingrecords((size_t)keyindex, value);
		}

		void printrecord(const size_t n){
			for (size_t mi = 0; mi < header.size(); mi++){
				std::printf("%s: %s\n", header[mi].c_str(), records[n][mi].c_str());
			}
			std::printf("\n");
		}
	};

class cMetaDataTable{

public:
	std::vector<std::string> header;
	std::vector<std::vector<std::string>> records;

	cMetaDataTable(){ }

	cMetaDataTable(const std::string metadatafile, size_t nskip){

		FILE* fp = fileopen(metadatafile, "r");

		std::string str;
		std::vector<std::string> tokens;

		//Read header line
		filegetline(fp, str);
		split(str, '\t', tokens);		
		header = tokens;
		size_t nfields = header.size();

		//Skip lines after header
		for (size_t i = 0; i < nskip; i++){
			filegetline(fp, str);
		}

		//Read remaining lines
		size_t k = 1 + nskip;
		while (filegetline(fp, str)){
			k++;
			tokens.resize(0);
			split(str, '\t', tokens);
			if (tokens.size() == nfields - 1){
				tokens.push_back("");
			}

			if (tokens.size() != nfields){
				std::printf("Error: On line %lu of file %s\n", k, metadatafile.c_str());
				std::printf("Error: The number of header items (%lu) does not match the number of data items (%lu)\n", nfields, tokens.size());
			}

			for (size_t i = 0; i < tokens.size(); i++){
				tokens[i] = trim(tokens[i]);
				tokens[i] = stripquotes(tokens[i]);
			}
			records.push_back(tokens);
		}
		fclose(fp);
	}

	cMetaDataTable(const cBlock& b, const std::string& id){
		std::vector<std::vector<std::string>> g = b.getblockleftright(id);
		for (size_t i = 0; i < g.size(); i++){
			addfield(g[i][0]);
			setfield(g[i][0], g[i][1]);
		}	
	}

	bool addfield(const std::string fname){
		header.push_back(fname);
		for (size_t i = 0; i < records.size(); i++){
			const std::string empty;
			records[i].push_back(empty);
		}
		return true;
	}

	bool setfield(const std::string fname, const std::string value, const size_t recindex){
		int k = findkeyindex(fname);
		if (k < 0)return false;
		records[recindex][(size_t)k] = value;
		return true;		
	}

	bool setfield(const std::string fname, const std::string value){
		int k = findkeyindex(fname);
		if (k < 0)return false;
		if (records.size() == 0){
			records.resize(1);
			records[0].resize(header.size());
		}

		for (size_t i = 0; i < records.size(); i++){
			records[i][(size_t)k] = value;
		}
		return true;
	}

	int findkeyindex(const std::string& fname){
		for (size_t i = 0; i < header.size(); i++){
			if (header[i] == fname){
				return (int)i;
			}
		}
		return -1;
	}

	std::vector<size_t> findmatchingrecords(const size_t keyindex, const int value){
		std::vector<size_t> indices;
		for (size_t i = 0; i < records.size(); i++){
			int n = atoi(records[i][keyindex].c_str());
			if (n == value){
				indices.push_back(i);
			}
		}
		return indices;
	}

	std::vector<size_t> findmatchingrecords(const std::string key, const int value)	{
		int keyindex = findkeyindex(key);
		if (keyindex < 0){
			std::vector<size_t> indices;
			return indices;
		}
		return findmatchingrecords((size_t)keyindex, value);
	}

	void printrecord(const size_t n){
		for (size_t mi = 0; mi < header.size(); mi++){
			std::printf("%s: %s\n", header[mi].c_str(), records[n][mi].c_str());
		}
		std::printf("\n");
	}
};

class cMetaDataRecord{

	public:
	std::vector<std::string> header;
	std::vector<std::string> values;

	cMetaDataRecord(){};

	cMetaDataRecord(const cBlock& b, const std::string& id){
		cMetaDataTable T(b, id);
		header = T.header;
		values = T.records[0];
	}

	int findkeyindex(const std::string& fname){
		for (size_t i = 0; i < header.size(); i++){
			if (header[i] == fname){
				return (int)i;
			}
		}
		return -1;
	}

	void print(){
		for (size_t mi = 0; mi < header.size(); mi++){
			std::printf("%s: %s\n", header[mi].c_str(), values[mi].c_str());
		}
		std::printf("\n");
	}
};

class cVariable{

	ILDataset* pDataset;

public:	
	std::string standardname;
	std::string longname;
	std::string aliaskey;
	std::string aliasvalue;
	std::vector<std::string> likelyfieldnames;
	std::vector<std::string> fieldnames;
	std::vector<std::string> requiredbythemes;
	std::string units;

	cVariable(cBlock b){
		standardname = b.getstringvalue("StandardName");
		longname = b.getstringvalue("LongName");		
		aliaskey = b.getstringvalue("SurveyInfoAlias");
		std::string s = b.getstringvalue("LikelyNames");
		likelyfieldnames = tokenize(s);
		units = b.getstringvalue("Units");
	}

	size_t setdataset(ILDataset* pDa){
		
		pDataset = pDa;
		fieldnames.resize(0);
		if (aliaskey.size() > 0 && aliaskey != cBlock::ud_string()){
			aliasvalue = pDataset->surveyinfofieldname(aliaskey);
			if (aliasvalue.size() > 0){
				fieldnames.push_back(aliasvalue);
			}
		}

		for (size_t i = 0; i < likelyfieldnames.size(); i++){
			std::string s = pDataset->fieldnamelike(likelyfieldnames[i]);
			if (s.size()>0){
				bool alreadyhavefield = false;
				for (size_t j = 0; j<fieldnames.size(); j++){
					if (s == fieldnames[j]){
						alreadyhavefield = true;
					}
					break;
				}
				if (alreadyhavefield==false){
					fieldnames.push_back(s);
				}
			}
		}
		return fieldnames.size();
	}

	bool hasaliaskey(){
		if (aliaskey == cBlock::ud_string()){
			return false;
		}
		return true;
	}

	bool aliasfieldexists(){
		if (hasaliaskey()){
			std::string s = pDataset->fieldnamelike(aliasvalue);
			if (s.size() > 0){
				return true;
			}
		}
		return false;
	}

	bool reportproblems(){

		bool status = true;
		if (hasaliaskey()){
			if (pDataset->hassurveyinfoid(aliaskey)){
				if (aliasfieldexists() == false){
					std::printf("Warning 3: %s\n", standardname.c_str());
					std::printf("\tAlias key %s is defined but field %s does not exist\n", aliaskey.c_str(), aliasvalue.c_str());
					status = false;
				}
			}			
		}		

		if (fieldnames.size() > 1){
			std::printf("Warning 4: %s\n", standardname.c_str());
			std::printf("\tMore than one likely field was found: ");
			for (size_t i = 0; i < fieldnames.size(); i++){
				std::printf("%s ", fieldnames[i].c_str());
			}
			std::printf("\n");
			status = false;
		}
		return status;

	}

	std::string likelyfield(){
		return fieldnames[0];
	}
};

class cTheme{

public:
	bool process;
	std::string Name;
	std::vector<std::vector<std::string>> required;

	cTheme(cBlock b){
		Name = b.getstringvalue("Name");
		process = b.getboolvalue("Process");
		if (process == false)return;

		std::vector<std::string> v = b.getblockstrings("RequiredVariables");
		for (size_t i = 0; i < v.size(); i++){
			required.push_back(tokenize(v[i]));
		}
	};
	
};

class cCRS{

public:
	std::string name = "";
	std::string epsg_code = "";
	double semi_major_axis = 0;
	double inverse_flattening = 0;
	bool valid = false;

	cCRS(std::string datum){
		if (datum == "GDA94"){
			name = datum;
			epsg_code = "EPSG:4283";
			semi_major_axis = 6378137.0;
			inverse_flattening = 298.257222101;
			valid = true;
		}
		else if (datum == "WGS84"){
			name = datum;
			epsg_code = "EPSG:4326";
			semi_major_axis = 6378137.0;
			inverse_flattening = 298.257223563;
			valid = true;
		}
		else if (datum == "AGD66"){
			name = datum;
			epsg_code = "EPSG:4326";
			semi_major_axis = 6378160.0;
			inverse_flattening = 298.25;
			valid = true;
		}
		else{
			name = datum;
		}
	};

};

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

		//std::vector<cBlock> bv = B.findblocks("Variable");
		//for (size_t i = 0; i < bv.size(); i++){
		//	cVariable v(bv[i]);
		//	Variables.push_back(v);
		//}

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

			//if (di % 10 != 0 && IDBName != "GSQP1029MAG"){
			//if (di%10 != 0) continue;
			
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
			
			NcFile ncFile(NCPath, NcFile::replace);
			add_lineindex(ncFile, D);
			add_groupbyline_variables(ncFile, D);			
			add_line_startend_points(ncFile, D);
			add_sample_variables(ncFile,D);
			add_global_metadata(ncFile, R);
			add_geospatial_metadata(ncFile, D);
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

	/*
	int findvariableindex(const std::string& fieldname)
	{
		_GSTITEM_
		for (size_t vi = 0; vi < Variables.size(); vi++){
			for (size_t k = 0; k < Variables[vi].likelyfieldnames.size(); k++){
				if (strcasecmp(Variables[vi].likelyfieldnames[k], fieldname) == 0){
					return (int)vi;
				}
			}
		}
		return -1;
	};*/

	void field_add_variable(ILField& F, NcFile& ncFile, const std::string& varname, NcVar& var, std::vector<NcDim> dims)
	{
		_GSTITEM_
		if (F.datatype().isbyte()){
			var = ncFile.addVar(varname, ncUbyte, dims);			
		}
		else if (F.datatype().isshort()){
			var = ncFile.addVar(varname, ncShort, dims);			
		}
		else if (F.datatype().isint()){
			var = ncFile.addVar(varname, ncInt, dims);			
		}
		else if (F.datatype().isfloat()){
			var = ncFile.addVar(varname, ncFloat, dims);			
		}
		else if (F.datatype().isdouble()){
			var = ncFile.addVar(varname, ncDouble, dims);			
		}
		else{
			std::string msg = "Error: Unknown Intrepid data type\n";
			message(flog,msg.c_str());
			throw(msg);
		}
	}

	void field_set_fillvalue(ILField& F, NcVar& var)
	{
		_GSTITEM_
		if (F.datatype().isbyte()){			
			var.putAtt("_FillValue", ncUbyte, IDatatype::bytenull());
		}
		else if (F.datatype().isshort()){			
			var.putAtt("_FillValue", ncShort, IDatatype::shortnull());
		}
		else if (F.datatype().isint()){			
			var.putAtt("_FillValue", ncInt, IDatatype::intnull());
		}
		else if (F.datatype().isfloat()){			
			var.putAtt("_FillValue", ncFloat, preferred_float_fillvalue());
		}
		else if (F.datatype().isdouble()){			
			var.putAtt("_FillValue", ncDouble, preferred_double_fillvalue());
		}
		else{
			std::string msg = "Error: Unknown Intrepid data type\n";
			message(flog, msg.c_str());
			throw(msg);
		}
	}

	bool add_lineindex(NcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		size_t nlines = D.nlines();
		NcDim  dim_line = ncFile.addDim("line", nlines);

		size_t nsamples = D.nsamples();
		NcDim  dim_sample = ncFile.addDim("sample", nsamples);

		NcVar vifirst = ncFile.addVar("index_first_sample", ncInt, dim_line);
		vifirst.putAtt("standard_name", "index_first_sample");
		vifirst.putAtt("description", "zero based index of the first sample in the line");
		vifirst.putAtt("units", "1");
		vifirst.putAtt("_FillValue", ncInt, IDatatype::intnull());

		NcVar vicount = ncFile.addVar("index_count", ncInt, dim_line);
		vicount.putAtt("standard_name", "index_count");
		vicount.putAtt("description", "number of samples in the line");
		vicount.putAtt("units", "1");
		vicount.putAtt("_FillValue", ncInt, IDatatype::intnull());

		std::vector<int> ifirst(nlines);
		std::vector<int> ilast(nlines);
		std::vector<int> icount(nlines);
		std::vector<int> isequence(nsamples);
		int startindex = 0;
		for (size_t li = 0; li < nlines; li++){
			ifirst[li] = startindex;
			icount[li] = (int)D.nsamplesinline(li);
			ilast[li] = ifirst[li] + icount[li] - 1;
			for (size_t k = (size_t)ifirst[li]; k <= (size_t)ilast[li]; k++){
				isequence[k] = (int)li;
			}
			startindex += icount[li];
		}
		vifirst.putVar(ifirst.data());		
		vicount.putVar(icount.data());		
		return true;
	}

	bool add_lineindex_unused(NcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		size_t nlines = D.nlines();

		size_t nsamples = D.nsamples();
		NcDim  dim_sample = ncFile.addDim("sample", nsamples);

		NcDim  dim_line = ncFile.addDim("line", nlines);
		NcDim  dim_line_index = ncFile.addDim("line_index", nlines+1);

		NcVar vifirst = ncFile.addVar("index_lines", ncInt, dim_line_index);
		vifirst.putAtt("standard_name", "index_lines");
		vifirst.putAtt("description", "zero based index of the first sample in the line, and then one past the end of the last line");
		vifirst.putAtt("units", "1");
		vifirst.putAtt("_FillValue", ncInt, IDatatype::intnull());

		std::vector<int> ifirst(nlines+1);
		std::vector<int> ilast(nlines);
		std::vector<int> icount(nlines);
		std::vector<int> isequence(nsamples);
		int startindex = 0;
		for (size_t li = 0; li < nlines; li++){
			ifirst[li] = startindex;
			icount[li] = (int)D.nsamplesinline(li);
			ilast[li]  = ifirst[li] + icount[li] - 1;
			for (size_t k = (size_t)ifirst[li]; k <= (size_t)ilast[li]; k++){
				isequence[k] = (int)li;
			}
			startindex += icount[li];
		}
		ifirst[nlines] = ilast[nlines-1] + 1;
		vifirst.putVar(ifirst.data());				
		return true;
	}

	bool add_line_startend_points(NcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		ILField& fx = D.getsurveyinfofield_ref("X");
		ILField& fy = D.getsurveyinfofield_ref("Y");

		size_t nlines = D.nlines();		
		NcDim  dim_line = ncFile.getDim("line");

		NcVar vx1 = ncFile.addVar("longitude_first", ncDouble, dim_line);
		vx1.putAtt("standard_name", "longitude_first");
		vx1.putAtt("description", "first non-null longitude coordinate in the line");
		vx1.putAtt("units", "degree_east");
		vx1.putAtt("_FillValue", ncDouble, preferred_double_fillvalue());

		NcVar vx2 = ncFile.addVar("longitude_last", ncDouble, dim_line);
		vx2.putAtt("standard_name", "longitude_last");
		vx2.putAtt("description", "last non-null longitude coordinate in the line");
		vx2.putAtt("units", "degree_east");
		vx2.putAtt("_FillValue", ncDouble, preferred_double_fillvalue());

		NcVar vy1 = ncFile.addVar("latitude_first", ncDouble, dim_line);
		vy1.putAtt("standard_name", "latitude_first");
		vy1.putAtt("description", "first non-null latitude coordinate in the line");
		vy1.putAtt("units", "degree_north");
		vy1.putAtt("_FillValue", ncDouble, preferred_double_fillvalue());

		NcVar vy2 = ncFile.addVar("latitude_last", ncDouble, dim_line);
		vy2.putAtt("standard_name", "latitude_last");
		vy2.putAtt("description", "last non-null latitude coordinate in the line");
		vy2.putAtt("units", "degree_north");
		vy2.putAtt("_FillValue", ncDouble, preferred_double_fillvalue());

		IDatatype dt = fx.datatype();

		std::vector<double> x1(nlines);
		std::vector<double> x2(nlines);
		std::vector<double> y1(nlines);
		std::vector<double> y2(nlines);		
		for (size_t li = 0; li < nlines; li++){
			ILSegment& sx = fx.Segments[li];
			ILSegment& sy = fy.Segments[li];
			sx.readbuffer();
			sy.readbuffer();
			size_t ns = sx.nsamples();			
			for(size_t k=0; k<ns; k++){
				x1[li] = sx.d(k);
				y1[li] = sy.d(k);
				if(dt.isnull(x1[li])==false && dt.isnull(y1[li])==false)break;
			}

			for (size_t k = ns-1; k!=0; k--){
				x2[li] = sx.d(k);
				y2[li] = sy.d(k);
				if (dt.isnull(x2[li]) == false && dt.isnull(y2[li]) == false)break;
			}					
		}

		vx1.putVar(x1.data());
		vx2.putVar(x2.data());
		vy1.putVar(y1.data());
		vy2.putVar(y2.data());
		return true;
	}

	bool add_groupbyline_variables(NcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;		
		size_t nlines   = D.nlines();				
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
			std::string units         = V.records[vi][i_units];						
						
			std::vector<NcDim> dims;
			dims.push_back(dim_line);
			if (F.nbands() > 1){
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			NcVar var;
			field_add_variable(F, ncFile, variable_name, var, dims);
			var.putAtt("standard_name", standard_name);			
			var.putAtt("original_database_name", F.Name);
			var.putAtt("units", units);
			field_set_fillvalue(F, var);
			
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
	
	bool add_sample_variables(NcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();		
		NcDim  dim_sample = ncFile.getDim("sample");
		
		int i_field         = V.findkeyindex("field_name");
		int i_convert       = V.findkeyindex("convert");
		int i_variable_name = V.findkeyindex("variable_name");
		int i_standard_name = V.findkeyindex("standard_name");
		int i_units         = V.findkeyindex("units");

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
			std::string units         = V.records[vi][i_units];						
			
			
			std::vector<NcDim> dims;
			dims.push_back(dim_sample);
			if (F.nbands() > 1){
				std::string dimname = "nbands_" + F.Name;
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			NcVar var;
			field_add_variable(F, ncFile, variable_name, var, dims);
			var.putAtt("standard_name", standard_name);
			var.putAtt("original_database_name", F.Name);		
			var.putAtt("units", units);									
			field_set_fillvalue(F, var);			
			if (DummyRun) continue;

			size_t startindex = 0;			
			for (size_t li = 0; li < nlines; li++){
				ILSegment& S = F.Segments[li];
				S.readbuffer();

				//Replace the nulls with more "sensible" values
				if (S.datatype().isfloat()) change_nulls_float(S, preferred_float_fillvalue());
				else if (S.datatype().isdouble()) change_nulls_double(S, preferred_double_fillvalue());

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

	bool add_global_metadata(NcFile& ncFile, const cMetaDataRecord& m){
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

	bool add_geospatial_metadata(NcFile& ncFile, ILDataset& D){
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
			cCRS c(datum);									
			
			//NcVar varcrs = ncFile.addVar("crs", ncInt);			
			std::vector<NcDim> dims;			
			int varid;	
			nc_def_var(ncFile.getId(),"crs",ncInt.getId(),0, NULL, &varid);
			NcVar varcrs(ncFile,varid);			

			varcrs.putAtt("grid_mapping_name", "latitude_longitude");
			varcrs.putAtt("epsg_code", c.epsg_code.c_str());
			varcrs.putAtt("semi_major_axis", ncDouble, c.semi_major_axis);
			varcrs.putAtt("inverse_flattening", ncDouble, c.inverse_flattening);
		}
		return true;
	}	

	static float preferred_float_fillvalue(){
		_GSTITEM_
		return -9999.0f;
	}

	static double preferred_double_fillvalue(){
		_GSTITEM_
		return -9999.0;
	}

	void change_nulls_float(ILSegment& S, const float& nullvalue){
		_GSTITEM_
		if (S.datatype().isfloat()){
			float* fp = (float*)S.pvoid();
			for (size_t k = 0; k < S.nstored(); k++){
				if (S.datatype().isnull(fp[k])){
					fp[k] = nullvalue;
				}
			}
		}
	}

	void change_nulls_double(ILSegment& S, const double& nullvalue){

		if (S.datatype().isdouble()){
			double* fp = (double*)S.pvoid();
			for (size_t k = 0; k < S.nstored(); k++){
				if (S.datatype().isnull(fp[k])){
					fp[k] = nullvalue;
				}
			}
		}

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


