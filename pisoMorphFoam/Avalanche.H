static scalarField getStencil(scalarField uArray, int start_index, int end_index)
{
    int uArraySize = uArray.size();
    int stencilSize = end_index - start_index;
    scalarField retval(stencilSize,0.0);
    
    // Case one - off the end..
    if (end_index > uArraySize)
    {
        int index = 0;
        for (int i=start_index; i < end_index; i++)
        {
            if (i < uArraySize)
            {
                retval[index] = uArray[i];
            }
            else
            {
                retval[index] = uArray[uArraySize - 1];
            }
            index +=1;
        }
    }
    else if (start_index < 0)
    {
        int index = 0;
        for (int i=start_index; i < end_index; i++)
        {
            if(i < 0)
            {
                retval[index] = uArray[0];
            }
            else
            {
                retval[index] = uArray[i];
            }
            index +=1;
        }
    }
    else
    {
        int index = 0;
        for (int i=start_index; i < end_index; i++)
        {
            retval[index] = uArray[i];
            index +=1;
        }
    }
    return retval;
}


// Returns the local slope in degrees
static scalarField getSlope(scalarField x, scalarField z)
{
	int n = x.size();
	scalarField slope(n);
    double dzdx = 0.0;
	forAll(x, index)
	{
        if (index == 0)
        {
            double dx = x[index + 1] - x[index];
            double dz = z[index + 1] - z[index];
            // Info << " x: " << x[index] << " z: " <<  z[index] << endl;
            if (dx > 0.0)
            {
                dzdx = dz/dx;
            }
            else
            {
                dzdx = 0.0;
            }
        }
        else if (index == n-1)
        {
            double dx = x[index] - x[index-1];
            double dz = z[index] - z[index-1];
            // Info << " x: " << x[index] << " z: " <<  z[index] << endl;
            if (dx > 0.0)
            {
                dzdx = dz/dx;
            }
            else
            {
                dzdx = 0.0;
            }
        }
		else
        {
            double dx = x[index + 1] - x[index - 1];
            double dz = z[index + 1] - z[index - 1];
            if (dx > 0.0)
            {
                dzdx = dz/dx;
            }
            else
            {
                dzdx = 0.0;
            }
            
        }
        double angle = Foam::atan(dzdx)*180.0/PI;
        slope[index] = angle;
	}
	return slope;
}

// This is the function that is solved for in the bisect method
static double func(double zTarget, scalarField xloc, scalarField zloc, double targetArea)
{
    double retval = 0.0;

    // cheap way of getting dx
    double dx = xloc[1] = xloc[0];
    int numPoints = 7;

    scalarField xtmp(numPoints,0.0);
    scalarField ztmp(numPoints,0.0);

    xtmp[0] = 0.0;
    xtmp[1] = 1.0*dx;
    xtmp[2] = 2.0*dx;
    xtmp[3] = dx*2.;
    xtmp[4] = dx;
    xtmp[5] = 0.0;
    xtmp[6] = xtmp[0];
        
    ztmp[0]=zloc[0];
    ztmp[1]=zTarget;
    ztmp[2]=zloc[2];
    ztmp[3]=0.0;
    ztmp[4]=0.0;
    ztmp[5]=0.0;
    ztmp[6]=zloc[0];

    int j = numPoints -1;

    for (int i=0; i < numPoints; i++)
    {
        retval += (xtmp[j] + xtmp[i]) * (ztmp[j] - ztmp[i]);
        j = i;
    }
    double area = retval/2.0;
    return area - targetArea;
}

// Returns the area within the polygon bounded by the arrays
// http://www.mathopenref.com/coordpolygonarea2.html
static double getAreaPolygon(scalarField xloc, scalarField zloc)
{
    double retval = 0.0;

    // cheap way of getting dx
    double dx = xloc[1] = xloc[0];
    int numPoints = 7;

    scalarField xtmp(numPoints,0.0);
    scalarField ztmp(numPoints,0.0);

    xtmp[0] = 0.0;
    xtmp[1] = 1.0*dx;
    xtmp[2] = 2.0*dx;
    xtmp[3] = 2.0*dx;
    xtmp[4] = 1.0*dx;
    xtmp[5] = 0.0;
    xtmp[6] = 0.0;
        
    ztmp[0]=zloc[0];
    ztmp[1]=zloc[1];
    ztmp[2]=zloc[2];
    ztmp[3]=0.0;
    ztmp[4]=0.0;
    ztmp[5]=0.0;
    ztmp[6]=zloc[0];

    int j = numPoints -1;

    for (int i=0; i < numPoints; i++)
    {
        retval += (xtmp[j] + xtmp[i]) * (ztmp[j] - ztmp[i]);
        j = i;
    }
    return retval/2.0;  
}


// https://www.geeksforgeeks.org/program-for-bisection-method/
static double bisect(scalarField xloc, scalarField zloc, double targetArea)
{
    double a = -1.0;
    double b = 100.;
    double c = a; 
    double TOL = 0.0001;

    int iteration = 0;
    
    while ((b-a) >= TOL && iteration < 100) 
    { 
        // Find middle point 
        c = (a+b)/2.; 
        double cError = func(c, xloc, zloc, targetArea);
        double aError = func(a, xloc, zloc, targetArea);
        // Check if middle point is root 
        if ( fabs(cError) < TOL)
	{
	    // Info << "Final bisect error " << fabs(cError) << " Iteration " << iteration << endl; 
            break; 
	}  
        // Decide the side to repeat the steps 
        else if ( cError*aError < 0) 
            b = c; 
        else
            a = c; 
        iteration += 1;
	// Info << "Bisect error " << cError << endl;
    } 
    return c;
}


// Iteratively adjusts the elevation so that the desired area is achieved.
static double adjustBed(scalarField xloc, scalarField zloc, double delArea)
{
    double currentArea = getAreaPolygon(xloc, zloc);
    double targetArea = currentArea + delArea;

    double zNew = bisect(xloc, zloc, targetArea);
    
    return zNew;
}


// Calculates the adjustment to the bed based on the angle os repose.
static double adjustToTarget(
        scalarField xloc, 
        scalarField zloc, 
        double angleOfRepose)
{
    double dx = 0.5*(xloc[2] - xloc[0]);

    double dzTarget = dx * Foam::tan(angleOfRepose*PI/180.);
    double dzCurrent = zloc[1]-zloc[2];
    
    return dzCurrent - dzTarget;
}

// Adjusts the profile using the Avalanche Scheme
static void avalancheProfile 
(
	vectorField& faceCentres,
        double thresholdAngle,
        double adjustmentAngle,
        double angleOfRepose
) 
{         
    int nSize = faceCentres.size(); 
    scalarField zArray(nSize);      // Current z elevations
    scalarField xArray(nSize);      // Current x locations
    scalarField zArrayNew(nSize);   // Array to hold changes to z
    

    int maxIterations = 100;    // Hard coding this 
    double factor = 0.65;         // Another hard coded value

    // Populate the coordinates into the scalarFields
    forAll(faceCentres, index)
    {
        xArray[index] = faceCentres[index].x();
        zArray[index] = faceCentres[index].z();
        zArrayNew[index] = faceCentres[index].z();
    }      
    

    scalarField slope = getSlope(xArray,zArray); 
    Info << "Max slope " << max(slope) << " Min slope " << min(slope) << endl; 

    Info << "-------------------------" << endl;
    Info << "Z Bed inital" << zArray << endl;
    Info << "-------------------------" << endl;
    
    // Check if the bed slope exceeds the threshold angle
    double minSlope = min(slope);
    int iterations = 0;
    if (minSlope < thresholdAngle* -1.0)
    {
        bool bedChanged = true;
        while (bedChanged == true && iterations < maxIterations)
        {
            bedChanged = false;            
            iterations += 1;
            Info << iterations << endl;
            // Main loop here
            forAll(xArray, index)
            {
                scalarField xloc = getStencil(xArray, index - 1, index + 2);
                scalarField zloc = getStencil(zArray, index - 1, index + 2);
                scalarField localSlope = getSlope(xloc, zloc);    
                
                if (min(localSlope) < angleOfRepose*-1.0)
                {

                    // Info << " Sub slope " << localSlope << endl;

                    double oldArea = getAreaPolygon(xloc, zloc);
                    // Calculate the adjustment to the z to match the angle of repose 
                    double dz = adjustToTarget(xloc, zloc, angleOfRepose);
                    
                    // This adjusts the rate at which we search
                    if (iterations < 75)
                    {   
                        zArrayNew[index] -= dz*factor;
                    }
                    else
                    {
                        zArrayNew[index] -= dz*0.65;
                    }
                    
                    // Get the new stencil 
                    scalarField zlocNew = getStencil(zArrayNew, index - 1, index + 2);
                    
                    // Calculate the new area (should be lower than previously)
                    double newArea = getAreaPolygon(xloc, zlocNew);
                    double deltaArea = oldArea - newArea;

                    // Info << "Old area: " << oldArea << " New Area: " << newArea << endl;
                    // Info << "deltaArea: " << deltaArea << endl;

                    // The mass moved from the current index needs to be moved
                    // downstream to the next node.
                    if (deltaArea > 0.)
                    {
                        // Get the i+1 stencil
                        scalarField xlocP1 = getStencil(xArray, index, index + 3);
                        scalarField zlocP1 = getStencil(zArray, index, index + 3);
                        
                        // iteratively search for the new z elevation to match the
                        // area to conserve the mass.
                        double zNew = adjustBed(xlocP1, zlocP1, deltaArea);
                        //  Note that I removed the periodic boundary adjustment
                        //  which pushes the upstream elevation to the downstream
                        zArrayNew[index + 1] = zNew; 
                    }
                    bedChanged = true;
                } // if min slope
            } // forAll   
            forAll(zArray, index)
            {
                zArray[index] = zArrayNew[index];
            }
        } // while
    } // end if on slope
    Info << "-------------------------" << endl;
    Info << "Z Bed Final" << zArray << endl;
    Info << "-------------------------" << endl;
    // return zArray;
} // end

