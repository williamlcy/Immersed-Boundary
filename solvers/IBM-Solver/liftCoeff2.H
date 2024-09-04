
// store the p interpolation at all IB points
std::vector<scalar> P_Interpolation(npt);
// Create an interpolation field based on the pressure field p
interpolationCell<scalar> pInterp(p);

// npt is the IB points number
for(label I=0; I<npt; I++)
{
    scalar P_Interpolation_Temp = -1000000000000.0; 
    vector aPoint = Lmarks[I];
    // find the cell include IB point
    label CellI = mesh.findCell(aPoint);
    // if the IB point is in the cell
    if (CellI != -1)
    {
        // cal the interpolation value at the IB point
        P_Interpolation_Temp = pInterp.interpolate
        (
            aPoint, 
            CellI
        );
    }
    
    // Pick max value from all processors
    reduce(P_Interpolation_Temp, maxOp<scalar>());
    P_Interpolation[I]=P_Interpolation_Temp;
}

// store the U interpolation at all IB points
std::vector<tensor> gradU_Interpolation(npt);
// Create an interpolation field based on the grad of U
gradU=fvc::grad(U);
interpolationCell<tensor> gradUInterp(gradU);

for(label I=0; I<npt; I++)
{
    tensor gradU_Interpolation_Temp(-1000000000000.0,-1000000000000.0,-1000000000000.0,
    -1000000000000.0,-1000000000000.0,-1000000000000.0,
    -1000000000000.0,-1000000000000.0,-1000000000000.0); 
    vector aPoint = Lmarks[I];
    label CellI = mesh.findCell(aPoint);
    if (CellI != -1)
    {
        gradU_Interpolation_Temp = gradUInterp.interpolate
        (
            aPoint, 
            CellI
        );
    }
    
    // Pick max value from all processors
    reduce(gradU_Interpolation_Temp, maxOp<tensor>());
    gradU_Interpolation[I]=gradU_Interpolation_Temp;
}

// pressure force and viscous force
vector pressureForce(0.0,0.0,0.0);
vector viscousForce(0.0,0.0,0.0);

for(label I=0; I<npt; I++)
{
    pressureForce += P_Interpolation[I]*Lmarks_normal[I]*Lmarks_Ds;
    //viscousForce += (gradU_Interpolation[I]&(Lmarks_tangential[I]))*(Lmarks_Ds*nu);
    auto oneV = (gradU_Interpolation[I]&(Lmarks_tangential[I]))*(Lmarks_Ds*nu);
    viscousForce += oneV.value();
}
Info << "Pressure Force by interpolationCellPoint: " << "F_x = " << pressureForce.x() << "  " << "F_y = " << pressureForce.y() << endl;
Info << "viscous Force by interpolationCellPoint: " << "F_x = " << viscousForce.x() << "  " << "F_y = " << viscousForce.y() << endl;
/*
if(myid==0)
{
getchar();
}
*/
MPI_Barrier(MPI_COMM_WORLD);
