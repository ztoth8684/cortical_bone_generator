// This is an extended version of the Measure command
// that calculates porosity of a CT scan

   setOption("Mean", true);
   addedporosity = 0;
   setSlice(1)
   
   for (i=1;i<nSlices+1; i++) {
   run("Measure");
   mean = getResult("Mean",nResults-1);
   porosity = 1- (mean/255);
   setResult('Porosity', nResults-1, porosity);
   run("Next Slice [>]");
   addedporosity = addedporosity + porosity;
   }
   
   totalporosity = addedporosity/nSlices;
   setResult('Tot.Porosity', 0, totalporosity);
   setResult('Tot.Porosity', nResults-1, totalporosity);
   updateResults();
