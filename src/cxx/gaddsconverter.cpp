
//Unused at the moment because we are not using GADDS data

class cGaddsConverter{

public:
	std::string LogFile;
	FILE* flog;
	std::vector<cTheme> Themes;
	std::vector<cVariable> Variables;
	cMetaDataTable M;
	cMetaDataTable A;
	cMetaDataTable J;
	std::string idbroot;
	std::string rmfromdodsurl;
	std::string ncoutdir;

	cGaddsConverter(const std::string& controlfile){
		cBlock B(controlfile);

		LogFile = B.getstringvalue("LogFile");		
		flog = fileopen(LogFile, "w");
		rootmessage(flog, "Program starting at %s\n", timestamp().c_str());
		B.write(flog);

		std::vector<cBlock> bv = B.findblocks("Theme");
		for (size_t i = 0; i < bv.size(); i++){
			cTheme t(bv[i]);
			if (t.process){
				Themes.push_back(t);
			}
		}

		bv = B.findblocks("Variable");
		for (size_t i = 0; i < bv.size(); i++){
			cVariable v(bv[i]);
			Variables.push_back(v);
		}

		M = cMetaDataTable(B, "Metadata");
		std::string argusfile = B.getstringvalue("ArgusMetaData");
		A = cMetaDataTable(argusfile, 1);

		std::string jetcatfile = B.getstringvalue("JetCatFile");
		J = cMetaDataTable(jetcatfile, 0);

		idbroot = B.getstringvalue("IntrepidDatabasesPath");
		rmfromdodsurl = B.getstringvalue("RemoveFromDODS_URL");
		ncoutdir = B.getstringvalue("NetCDFOutputDir");		
	};

	~cGaddsConverter(){
		rootmessage("Finished\n");
		rootmessage(flog, "Finished at %s\n", timestamp().c_str());
		fclose(flog);
	};

	bool process(){
		
		bool status = checkthemes();
		if (status == false)return status;

		int j_dodsurl = J.findkeyindex("DODS_URL");
		int j_datatype = J.findkeyindex("DATATYPE");
		int j_theme = J.findkeyindex("THEME");
		int j_surveyid = J.findkeyindex("Surveyid");
		for (size_t ji = 0; ji < J.records.size(); ji++){
			//printf("Jetcat Record %lu\n", ji);
			size_t jetcatrecord = ji;
			std::string upath = J.records[ji][(size_t)j_dodsurl];
			std::string dtype = J.records[ji][(size_t)j_datatype];
			std::string theme = J.records[ji][(size_t)j_theme];
			std::string surveyid = J.records[ji][(size_t)j_surveyid];

			if (theme != "GRAVITY" && theme != "RADIOMETRICS" && 			theme != "ELEVATION" && theme != "MAGNETICS"){
				message(flog, "Warning unknown theme <%s> in record %lu url=%s\n", theme.c_str(), ji, upath.c_str());
				continue;
			}

			if (theme == "GRAVITY"){
				//message(flog,"Warning skipping GRAVITY dataset 				in record %lu url=%s\n", di, upath.c_str());
				continue;
			}

			if (dtype != "VECTOR"){
				//message(flog,"Warning skipping non VECTOR dataset 				in record %lu url=%s\n", ji, upath.c_str());
				continue;
			}

			int surveynumber = atoi(surveyid.c_str());
			std::vector<size_t> arecords = A.findmatchingrecords("PROJECT", surveynumber);
			if (arecords.size() == 0){
				message(flog,"Warning could not find matching PROJECT for Jetcat record %lu url=%s\n", ji, upath.c_str());
				continue;
			}
			else if (arecords.size() > 1){
				message(flog, "Warning found more than 1 matching PROJECT for Jetcat record %lu url=%s\n", ji, upath.c_str());
				continue;
			}
			size_t argusrecord = arecords[0];

			cMetaDataTable T = M;
			for (size_t i = 0; i < T.header.size(); i++){
				std::string& s = T.records[0][i];
				if (std::strncmp(s.c_str(), "Argus.", 6) == 0){
					std::string argusfield = s.substr(6, s.size() - 6);
					int ki = A.findkeyindex(argusfield);
					s = A.records[argusrecord][ki];
				}

				if (std::strncmp(s.c_str(), "Jetcat.", 7) == 0){
					std::string jetcatfield = s.substr(7, s.size() - 7);
					int ki = J.findkeyindex(jetcatfield);
					s = J.records[jetcatrecord][ki];
				}

				int ki = T.findkeyindex("date_created");
				if (ki >= 0){
					//std::string datestr = timestring("%Y%m%d");
					std::string datestr = timestring("%Y-%m-%d %H:%M:%S");
					T.records[0][ki] = datestr;
				}
			}
			//T.printrecord(0);

			size_t pos = upath.find(rmfromdodsurl) + rmfromdodsurl.size();
			std::string dspath = idbroot + upath.substr(pos, upath.size() - pos);
			fixseparator(dspath);

			bool datasetexists = exists(dspath);
			if (datasetexists == false){
				message(flog,"Warning: dataset for record %lu does not exist: %s\n", ji, dspath.c_str());
				continue;
			}

			sFilePathParts fpp = getfilepathparts(dspath);
			pos = dspath.find("..DIR");
			std::string dbpath = dspath.substr(0, pos);
			std::string dbname = extractfilename(dspath);

			std::string ncpath = ncoutdir + dbname + ".nc";
			std::vector<cMetaDataRecord> metadata(2);
			
			//std::string s = "\\\\nas\\cdsm\\application\\GADDS\\wa\\GSWA_P1261MAG";
			//if (s != dbpath)continue;
			message(flog,"%lu: Converting database %s\n", ji, dbpath.c_str());			
			checkdatabase(dbpath, theme);

			

			//double t1   = gettime();			
			//bool status = fieldstats(dbpath);
			//double t2   = gettime();
			//std::printf("time=%lf\n\n", t2 - t1);

			//A.printrecord(argusrecord);
			//std::printf("\n");
			//J.printrecord(jetcatrecord);
			//std::printf("\n");
			//intrepid2netcdf_flat(dbpath, ncpath, metadata);
		}
		return true;
	}

	int themeindex(const std::string themename){
		for (size_t i = 0; i < Themes.size(); i++){
			if (themename == Themes[i].Name){
				return (int)i;
			}
		}
		return -1;
	}

	int variableindex(const std::string variablename){
		for (size_t i = 0; i < Variables.size(); i++){
			if (variablename == Variables[i].standardname){
				return (int)i;
			}
		}
		return -1;
	}

	bool checkthemes(){

		bool status = true;
		for (size_t ti = 0; ti < Themes.size(); ti++){
			for (size_t ri = 0; ri < Themes[ti].required.size(); ri++){
				std::vector<std::string> rvars = Themes[ti].required[ri];
				for (size_t vi = 0; vi < rvars.size(); vi++){
					int k = variableindex(rvars[vi]);
					if (k >= 0)continue;
					message(flog, "Error: Required variable %s in theme %s is used but was not defined\n",rvars[vi].c_str(),Themes[ti].Name.c_str());
					status = false;
				}
			}
		}
		return status;
	}

	bool checkdatabase(const std::string& dbpath, const std::string& themename){		
		
		int ti = themeindex(themename);
		if (ti < 0)return false;

		class ILDataset D(dbpath);
		if (D.valid == false){
			message(flog, "Error: invalid database %s", dbpath.c_str());
			return false;
		}

		cTheme& theme = Themes[ti];
		for (size_t ri = 0; ri < theme.required.size(); ri++){						
			
			int nfields = 0;
			std::vector<std::string> vars = theme.required[ri];
			for (size_t i = 0; i < vars.size(); i++){
				int vi = variableindex(vars[i]);
				if (vi < 0)continue;
				cVariable& var = Variables[vi];
				size_t nf = var.setdataset(&D);
				nfields += nf;
				bool hstatus = var.reportproblems();
				//if (hstatus == false){
				//	D.printinfo();
				//	std::printf("Error not find a field for %s\n", var.standardname.c_str());
				//	prompttocontinue();
				//}				
			}

			if (nfields == 0){
				message(flog, D.infostring().c_str());				
				message(flog,"Error: Could not find a field for any of the variables: ");
				for (size_t k = 0; k < vars.size(); k++){
					message(flog, "%s ", vars[k].c_str());
				}
				message(flog, "\n");				
			}
		}

	}
};
