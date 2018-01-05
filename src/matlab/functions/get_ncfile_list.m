function F = get_ncfile_list(catalog);

baseurl = 'http://dapds00.nci.org.au/thredds/';

url = [baseurl catalog];
s = urlread(url);

a = strfind(s,'dataset=');
k=1;
for i=2:1:length(a)
    b = strfind(s(a(i):end),'.nc');
    k1 = a(i)+8;
    k2 = k1+b(1)-7;
    name = s(k1:k2);
    F(k).ncurl = [baseurl 'dodsC/' name];    
    %disp(F(k).ncurl)
    k=k+1;
end