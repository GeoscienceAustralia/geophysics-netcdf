
int convertaseggdf2netcdf(const std::string& dfnfile){

	double t1, t2;

	sFilePathParts fpp = getfilepathparts(dfnfile);
	std::string datfilename = fpp.directory + fpp.prefix + ".dat";
	std::string ncfilename = fpp.directory + fpp.prefix + ".nc";


	cGeophysicsNcFile ncFile(ncfilename, NcFile::replace);

	t1 = gettime();
	size_t npoints = countlines(datfilename);
	t2 = gettime();
	printf("Number of points is %d time was %lf\n", npoints, t2 - t1);

	cAsciiColumnFile A;
	printf("Opening data file\n");
	A.openfile(datfilename);
	printf("Started parsing ASEGGDF2 header\n");
	A.parse_aseggdf2_header(dfnfile);
	printf("Finished parsing ASEGGDF2 header\n");

	size_t nfields = A.fields.size();
	//size_t fi_survey = A.fieldindexbyname("project_ga");
	//size_t fi_date = A.fieldindexbyname("date");
	//size_t fi_flight = A.fieldindexbyname("flight");
	size_t fi_line = A.fieldindexbyname("line");
	size_t fi_x = A.fieldindexbyname("East_MGA51");
	size_t fi_y = A.fieldindexbyname("North_MGA51");

	t1 = gettime();
	std::vector<std::vector<int>>    intfields;
	std::vector<std::vector<double>> dblfields;

	NcDim  dim_points = ncFile.addDim("points", npoints);
	//NcType type_integer = ncFile.addVlenType("vlen_integer", ncInt);
	//NcType type_double  = ncFile.addVlenType("vlen_double", ncDouble);

	std::vector<NcVar> var_field(nfields);
	for (size_t fi = 0; fi < nfields; fi++){
		//printf("Pre processing field %d\n", (int)fi);
		std::string fieldname = A.fields[fi].name;
		size_t nbands = A.fields[fi].nbands;

		vector<NcDim> dims;
		dims.push_back(dim_points);
		if (nbands > 1){
			std::string dimname = "nbands_" + fieldname;
			NcDim dimband = ncFile.addDim(dimname, nbands);
			dims.push_back(dimband);
		}

		if (A.fields[fi].isinteger()){
			var_field[fi] = ncFile.addVar(A.fields[fi].name, ncInt, dims);
			var_field[fi].putAtt("nullvalue", ncInt, A.fields[fi].nullvalue);
		}
		else{
			var_field[fi] = ncFile.addVar(A.fields[fi].name, ncDouble, dims);
			var_field[fi].putAtt("nullvalue", ncDouble, A.fields[fi].nullvalue);
		}

		var_field[fi].putAtt("long_name", fieldname);
		var_field[fi].putAtt("standard_name", fieldname);
		var_field[fi].putAtt("units", fieldname);

	}

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	size_t startindex = 0;
	size_t lineindex = 0;
	while (size_t nsamples = A.readnextgroup(fi_line, intfields, dblfields)){
		std::printf("Processing line %d\n", lineindex + 1);
		for (size_t fi = 0; fi < nfields; fi++){
			//printf("Writing field %d\n", fi);
			std::string fieldname = A.fields[fi].name;
			size_t nbands = A.fields[fi].nbands;

			std::vector<size_t> startp(2);
			std::vector<size_t> countp(2);
			startp[0] = startindex;
			countp[0] = nsamples;

			for (size_t bi = 0; bi < nbands; bi++){
				startp[1] = bi;
				countp[1] = 1;

				if (A.fields[fi].isinteger()){
					std::vector<int> data(nsamples);
					for (size_t si = 0; si < nsamples; si++){
						data[si] = intfields[fi][si*nbands + bi];
					}
					var_field[fi].putVar(startp, countp, data.data());
				}
				else{
					std::vector<double> data(nsamples);
					for (size_t si = 0; si < nsamples; si++){
						data[si] = dblfields[fi][si*nbands + bi];
					}
					var_field[fi].putVar(startp, countp, data.data());
				}
			}

			vector<double> east = dblfields[fi_x];
			vector<double> north = dblfields[fi_y];
			vector<double> longitude;
			vector<double> latitude;
			bool status = projection(east, north, longitude, latitude);

			auto xmme = std::minmax_element(longitude.begin(), longitude.end());
			auto ymme = std::minmax_element(latitude.begin(), latitude.end());
			if (*xmme.first  < xmin)xmin = *xmme.first;
			if (*xmme.second > xmax)xmax = *xmme.second;
			if (*ymme.first  < ymin)ymin = *ymme.first;
			if (*ymme.second > ymax)ymax = *ymme.second;
		}
		startindex += nsamples;
		lineindex++;
	}
	A.closefile();
	t2 = gettime();
	printf("Reading data elapsed time = %.2lf\n", t2 - t1);

	NcDim dimcontainer = ncFile.addDim("container", 1);
	NcVar varcrs = ncFile.addVar("crs", ncInt, dimcontainer);
	varcrs.putAtt("grid_mapping_name", "latitude_longitude");
	varcrs.putAtt("epsg_code", "EPSG:4283");
	varcrs.putAtt("semi_major_axis", ncDouble, 6378137.0);
	varcrs.putAtt("inverse_flattening", ncDouble, 298.257222101);

	ncFile.putAtt("Conventions", "CF - 1.6"); //......................................... REQUIRED    - Always try to use latest value. (CF)
	ncFile.putAtt("Metadata_Conventions", "Unidata Dataset Discovery v1.0"); //........ REQUIRED    - Do not change. (ACDD)
	ncFile.putAtt("featureType", "trajectory"); //..................................... REQUIRED - CF attribute for identifying the featureType.
	ncFile.putAtt("cdm_data_type", "Trajectory"); //................................... REQUIRED (ACDD)
	ncFile.putAtt("nodc_template_version", "NODC_NetCDF_Trajectory_Template_v1.1"); //....... REQUIRED (NODC)
	ncFile.putAtt("standard_name_vocabulary", "NetCDF Climate and Forecast(CF) Metadata Convention Standard Name Table X"); //........ REQUIRED    - If using CF standard name attribute for variables. "X" denotes the table number  (ACDD)

	ncFile.putAtt("title", "Southern Thomson Airborne Electromagnetic Survey 2015 - Layered Earth Inversions");
	ncFile.putAtt("summary", "The dataset is the point located conductivity model derived from the GALEISBSTDEM inversion program");
	ncFile.putAtt("keywords", "airborne, geophysics, electromagnetic, AEM, inversion, layered earth, conductivity model, Southern Thomson Orogen");
	ncFile.putAtt("id", "Geocat#1257");
	ncFile.putAtt("naming authority", "Geoscience Australia");
	ncFile.putAtt("institution", "Geoscience Australia");
	ncFile.putAtt("project", "Geoscience Australia Precompetitive Data Acquition");

	ncFile.putAtt("geospatial_lon_min", ncDouble, xmin);
	ncFile.putAtt("geospatial_lon_max", ncDouble, xmax);
	ncFile.putAtt("geospatial_lon_units", "degrees_east");
	ncFile.putAtt("geospatial_lon_resolution", "point");

	ncFile.putAtt("geospatial_lat_min", ncDouble, ymin);
	ncFile.putAtt("geospatial_lat_max", ncDouble, ymax);
	ncFile.putAtt("geospatial_lat_units", "degrees_north");
	ncFile.putAtt("geospatial_lat_resolution", "point");

	ncFile.putAtt("geospatial_vertical_min", ncDouble, 0.0);
	ncFile.putAtt("geospatial_vertical_max", ncDouble, 0.0);
	ncFile.putAtt("geospatial_vertical_units", "m");
	ncFile.putAtt("geospatial_vertical_resolution", "point");
	ncFile.putAtt("geospatial_vertical_positive", "up");
	return 0;
}

int convertaseggdf2netcdf_vlen(const std::string& dfnfile){

	sFilePathParts fpp = getfilepathparts(dfnfile);
	std::string datfilename = fpp.directory + fpp.prefix + ".dat";
	std::string ncfilename = fpp.directory + fpp.prefix + ".vlen.h5";

	NcFile ncFile(ncfilename, NcFile::replace);

	cAsciiColumnFile A;
	printf("Opening data file\n");
	A.openfile(datfilename);
	printf("Started parsing ASEGGDF2 header\n");
	A.parse_aseggdf2_header(dfnfile);
	printf("Finished parsing ASEGGDF2 header\n");

	size_t nfields = A.fields.size();
	//size_t fi_survey = A.fieldindexbyname("project_ga");
	//size_t fi_date = A.fieldindexbyname("date");
	//size_t fi_flight = A.fieldindexbyname("flight");
	size_t fi_line = A.fieldindexbyname("line");
	size_t fi_x = A.fieldindexbyname("emloop_easting");
	size_t fi_y = A.fieldindexbyname("emloop_northing");

	double t1, t2;
	t1 = gettime();
	std::vector<std::vector<int>>    intfields;
	std::vector<std::vector<double>> dblfields;

	size_t nlines = 101;
	NcDim  dim_lines = ncFile.addDim("trajectory", nlines);
	NcType type_vlen_integer = ncFile.addVlenType("vlen_integer", ncInt);
	NcType type_vlen_double = ncFile.addVlenType("vlen_double", ncDouble);

	std::vector<NcVar> var_field(nfields);
	for (size_t fi = 0; fi < nfields; fi++){
		//printf("Pre processing field %d\n", (int)fi);
		std::string fieldname = A.fields[fi].name;
		size_t nbands = A.fields[fi].nbands;

		vector<NcDim> dims;
		dims.push_back(dim_lines);
		if (nbands > 1){
			std::string dimname = "nbands_" + fieldname;
			NcDim dimband = ncFile.addDim(dimname, nbands);
			dims.push_back(dimband);
		}

		if (A.fields[fi].isinteger()){
			var_field[fi] = ncFile.addVar(A.fields[fi].name, type_vlen_integer, dims);
			var_field[fi].putAtt("nullvalue", ncInt, A.fields[fi].nullvalue);
		}
		else{
			var_field[fi] = ncFile.addVar(A.fields[fi].name, type_vlen_double, dims);
			var_field[fi].putAtt("nullvalue", ncDouble, A.fields[fi].nullvalue);
		}

		var_field[fi].putAtt("long_name", fieldname);
		var_field[fi].putAtt("standard_name", fieldname);
		var_field[fi].putAtt("units", fieldname);

	}

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	int li = 0;
	while (size_t nsamples = A.readnextgroup(fi_line, intfields, dblfields)){
		printf("Processing line %d\n", li + 1);
		if (li >= nlines)break;
		for (size_t fi = 0; fi < nfields; fi++){
			//printf("Writing field %d\n", fi);
			std::string fieldname = A.fields[fi].name;
			size_t nbands = A.fields[fi].nbands;

			std::vector<size_t> startp(2);
			std::vector<size_t> countp(2);
			startp[0] = li;
			countp[0] = 1;

			for (size_t bi = 0; bi < nbands; bi++){
				startp[1] = bi;
				countp[1] = 1;

				nc_vlen_t vldata;
				if (A.fields[fi].isinteger()){
					std::vector<int> data(nsamples);
					for (size_t si = 0; si < nsamples; si++){
						data[si] = intfields[fi][si*nbands + bi];
					}
					vldata.len = nsamples;
					vldata.p = data.data();
					var_field[fi].putVar(startp, countp, &vldata);
				}
				else{
					std::vector<double> data(nsamples);
					for (size_t si = 0; si < nsamples; si++){
						data[si] = dblfields[fi][si*nbands + bi];
					}
					vldata.len = nsamples;
					vldata.p = data.data();
					var_field[fi].putVar(startp, countp, &vldata);
				}
			}

			vector<double> east = dblfields[fi_x];
			vector<double> north = dblfields[fi_y];
			vector<double> longitude;
			vector<double> latitude;
			bool status = projection(east, north, longitude, latitude);

			auto xmme = std::minmax_element(longitude.begin(), longitude.end());
			auto ymme = std::minmax_element(latitude.begin(), latitude.end());
			if (*xmme.first  < xmin)xmin = *xmme.first;
			if (*xmme.second > xmax)xmax = *xmme.second;
			if (*ymme.first  < ymin)ymin = *ymme.first;
			if (*ymme.second > ymax)ymax = *ymme.second;
		}
		li++;
	}
	A.closefile();
	t2 = gettime();
	printf("Reading data elapsed time = %.2lf\n", t2 - t1);

	NcDim dimcontainer = ncFile.addDim("container", 1);
	NcVar varcrs = ncFile.addVar("crs", ncInt, dimcontainer);
	varcrs.putAtt("grid_mapping_name", "latitude_longitude");
	varcrs.putAtt("epsg_code", "EPSG:4283");
	varcrs.putAtt("semi_major_axis", ncDouble, 6378137.0);
	varcrs.putAtt("inverse_flattening", ncDouble, 298.257222101);

	ncFile.putAtt("Conventions", "CF - 1.6"); //......................................... REQUIRED    - Always try to use latest value. (CF)
	ncFile.putAtt("Metadata_Conventions", "Unidata Dataset Discovery v1.0"); //........ REQUIRED    - Do not change. (ACDD)
	ncFile.putAtt("featureType", "trajectory"); //..................................... REQUIRED - CF attribute for identifying the featureType.
	ncFile.putAtt("cdm_data_type", "Trajectory"); //................................... REQUIRED (ACDD)
	ncFile.putAtt("nodc_template_version", "NODC_NetCDF_Trajectory_Template_v1.1"); //....... REQUIRED (NODC)
	ncFile.putAtt("standard_name_vocabulary", "NetCDF Climate and Forecast(CF) Metadata Convention Standard Name Table X"); //........ REQUIRED    - If using CF standard name attribute for variables. "X" denotes the table number  (ACDD)

	ncFile.putAtt("title", "Southern Thomson Airborne Electromagnetic Survey 2015 - Layered Earth Inversions");
	ncFile.putAtt("summary", "The dataset is the point located conductivity model derived from the GALEISBSTDEM inversion program");
	ncFile.putAtt("keywords", "airborne, geophysics, electromagnetic, AEM, inversion, layered earth, conductivity model, Southern Thomson Orogen");
	ncFile.putAtt("id", "Geocat#1257");
	ncFile.putAtt("naming authority", "Geoscience Australia");
	ncFile.putAtt("institution", "Geoscience Australia");
	ncFile.putAtt("project", "Geoscience Australia Precompetitive Data Acquition");

	ncFile.putAtt("geospatial_lon_min", ncDouble, xmin);
	ncFile.putAtt("geospatial_lon_max", ncDouble, xmax);
	ncFile.putAtt("geospatial_lon_units", "degrees_east");
	ncFile.putAtt("geospatial_lon_resolution", "point");

	ncFile.putAtt("geospatial_lat_min", ncDouble, ymin);
	ncFile.putAtt("geospatial_lat_max", ncDouble, ymax);
	ncFile.putAtt("geospatial_lat_units", "degrees_north");
	ncFile.putAtt("geospatial_lat_resolution", "point");

	ncFile.putAtt("geospatial_vertical_min", ncDouble, 0.0);
	ncFile.putAtt("geospatial_vertical_max", ncDouble, 0.0);
	ncFile.putAtt("geospatial_vertical_units", "m");
	ncFile.putAtt("geospatial_vertical_resolution", "point");
	ncFile.putAtt("geospatial_vertical_positive", "up");
	return 0;
}

int convertaseggdf2netcdf_grouped(const std::string& dfnfile){

	sFilePathParts fpp = getfilepathparts(dfnfile);
	std::string datfilename = fpp.directory + fpp.prefix + ".dat";
	std::string ncfilename = fpp.directory + fpp.prefix + ".grouped.h5";

	NcFile ncFile(ncfilename, NcFile::replace);

	cAsciiColumnFile A;
	printf("Opening data file\n");
	A.openfile(datfilename);
	printf("Started parsing ASEGGDF2 header\n");
	A.parse_aseggdf2_header(dfnfile);
	printf("Finished parsing ASEGGDF2 header\n");
	size_t nfields = A.fields.size();

	size_t fi_survey = A.fieldindexbyname("project_ga");
	size_t fi_date = A.fieldindexbyname("date");
	size_t fi_flight = A.fieldindexbyname("flight");
	size_t fi_line = A.fieldindexbyname("line");
	size_t fi_x = A.fieldindexbyname("emloop_easting");
	size_t fi_y = A.fieldindexbyname("emloop_northing");

	double t1, t2;
	t1 = gettime();
	std::vector<std::vector<int>>    intfields;
	std::vector<std::vector<double>> dblfields;


	std::vector<NcDim> banddim(nfields);
	for (size_t fi = 0; fi < nfields; fi++){
		//printf("Pre processing field %d\n", (int)fi);		
		if (A.fields[fi].nbands > 1){
			std::string dimname = "nbands_" + A.fields[fi].name;
			banddim[fi] = ncFile.addDim(dimname, A.fields[fi].nbands);
		}
	}

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	int gi = 0;
	while (size_t nsamples = A.readnextgroup(fi_line, intfields, dblfields)){
		gi++;
		if (gi > 3)break;

		printf("Processing group %d\n", gi);
		std::string groupname = strprint("Line:%d", gi);
		NcGroup group = ncFile.addGroup(groupname);
		NcDim DimSamples = group.addDim("nsamples", nsamples);
		group.putAtt("SurveyNumber", ncInt, intfields[fi_survey][0]);
		group.putAtt("Date", ncInt, intfields[fi_date][0]);
		group.putAtt("FlightNumber", ncInt, intfields[fi_flight][0]);
		group.putAtt("LineNumber", ncInt, intfields[fi_line][0]);

		for (size_t fi = 0; fi < nfields; fi++){

			if (fi == fi_survey)continue;
			else if (fi == fi_date)continue;
			else if (fi == fi_flight)continue;
			else if (fi == fi_line)continue;

			//printf("Processing field %d\n", fi);
			std::string fieldname = A.fields[fi].name;
			size_t nbands = A.fields[fi].nbands;

			std::vector<NcDim> vardims;
			vardims.push_back(DimSamples);
			if (nbands > 1)vardims.push_back(banddim[fi]);

			NcVar field;
			if (A.fields[fi].isinteger()){
				field = group.addVar(fieldname, ncInt, vardims);
				field.putVar(intfields[fi].data());
				field.putAtt("nullvalue", ncInt, (int)A.fields[fi].nullvalue);
			}
			else if (A.fields[fi].isreal()){
				field = group.addVar(fieldname, ncDouble, vardims);
				field.putVar(dblfields[fi].data());
				field.putAtt("nullvalue", ncDouble, (int)A.fields[fi].nullvalue);
			}
			field.putAtt("long_name", fieldname);
			field.putAtt("standard_name", fieldname);
			field.putAtt("units", A.fields[fi].units);


			auto xmme = std::minmax_element(dblfields[fi_x].begin(), dblfields[fi_x].end());
			auto ymme = std::minmax_element(dblfields[fi_y].begin(), dblfields[fi_y].end());
			if (*xmme.first  < xmin)xmin = *xmme.first;
			if (*xmme.second > xmax)xmax = *xmme.second;
			if (*ymme.first  < ymin)ymin = *ymme.first;
			if (*ymme.second > ymax)ymax = *ymme.second;

		}
	}
	A.closefile();
	t2 = gettime();
	printf("Reading data elapsed time = %.2lf\n", t2 - t1);

	ncFile.putAtt("title", "Southern Thomson Airborne Electromagnetic Survey 2015 - Layered Earth Inversions");
	ncFile.putAtt("summary", "The dataset is the point located conductivity model derived from the GALEISBSTDEM inversion program");
	ncFile.putAtt("keywords", "airborne, geophysics, electromagnetic, AEM, inversion, layered earth, conductivity model, Southern Thomson Orogen");
	ncFile.putAtt("conventions", "CF, ACDD");
	ncFile.putAtt("id", "Geocat#1257");
	ncFile.putAtt("naming authority", "Geoscience Australia");
	ncFile.putAtt("institution", "Geoscience Australia");
	ncFile.putAtt("project", "Geoscience Australia Precompetitive Data Acquition");

	ncFile.putAtt("lon_min", ncDouble, xmin);
	ncFile.putAtt("lon_max", ncDouble, xmax);
	ncFile.putAtt("lat_min", ncDouble, ymin);
	ncFile.putAtt("lat_max", ncDouble, ymax);	
	return 0;
}

int testvlen()
{
	std::string ncfilename = "Z:\\code\\test_data\\aseggdf2netcdf\\test.h5";
	NcFile ncFile(ncfilename, NcFile::replace);

	static int nlines = 3;
	int nbands = 5;
	NcDim dimlines = ncFile.addDim("nlines", nlines);
	NcDim dimbands = ncFile.addDim("nbands", nbands);
	NcVlenType vlentype = ncFile.addVlenType("vlenType_1", ncDouble);

	vector<NcDim> dims(2);
	dims[0] = dimlines;
	dims[1] = dimbands;
	NcVar field = ncFile.addVar("field", vlentype, dims);

	for (size_t li = 0; li < nlines; li++){
		int nsamples = li + 4;
		nc_vlen_t vldata;
		for (size_t bi = 0; bi < nbands; bi++){
			std::vector<size_t> start(2);
			std::vector<size_t> count(2);
			start[0] = li;
			start[1] = bi;
			count[0] = 1;
			count[1] = 1;

			vector<double> d(nsamples);
			for (size_t si = 0; si < nsamples; si++){
				d[si] = li * 100 + bi * 10 + si;
			}
			vldata.len = d.size();
			vldata.p = d.data();
			field.putVar(start, count, &vldata);
		}
	}
	return 0;
}

