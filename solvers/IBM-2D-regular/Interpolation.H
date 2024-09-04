
std::vector<vector> U_Interpolation;
//std::vector<scalar> P_Interpolation;

/*
for(int I=0; I<Mag; I++)
{
    int NN = meshSupport[0][I];
    Info << "No." << I+1 << ": " << NN << endl;
    Info << "Velocity: " << "Ux = " << U[NN].component(0) << "  " << "Uy = " << U[NN].component(1) << endl;
}

dictionary interpolationDict = mesh.schemesDict().subDict("interpolationSchemes");
Info << "interpolationSchemes: " << interpolationDict << endl;
autoPtr<interpolation<scalar> > P_int = interpolation<scalar>::New(interpolationDict, p);
*/
std::vector<std::vector<double>> velLocal_Collect[3];
{

    // get the velocity of the local processor
    std::vector<double> velLocal[3];
    velLocal[0].resize(local_cellnum);
    velLocal[1].resize(local_cellnum);
    velLocal[2].resize(local_cellnum);
    for(int j=0;j<local_cellnum;j++)
    {
        velLocal[0][j]=U[j].component(0);
        velLocal[1][j]=U[j].component(1);
        velLocal[2][j]=U[j].component(2);
    }

    //Local velocity data is collected into a global container
    velLocal_Collect[0].resize(mpi_num);
    velLocal_Collect[1].resize(mpi_num);
    velLocal_Collect[2].resize(mpi_num);
    for(int k=0;k<mpi_num;k++)
    {
            velLocal_Collect[0][k].resize(cellNum[k]);
        velLocal_Collect[1][k].resize(cellNum[k]);
            velLocal_Collect[2][k].resize(cellNum[k]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(int ip=0;ip<3;ip++)
    {

        //send the velocity data of the local processor to the root processor
        if(myid!=0)
        {
            MPI_Send(velLocal[ip].data(), local_cellnum, MPI_DOUBLE, 0, 555, MPI_COMM_WORLD);
        }
        else
        {
            velLocal_Collect[ip][0]=velLocal[ip];
            for(int k=1;k<mpi_num;k++)
            {
             MPI_Recv(velLocal_Collect[ip][k].data(), cellNum[k], MPI_DOUBLE, k, 555, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //send the velocity data of the root processor to all the processors
        for(int kk=0;kk<mpi_num;kk++)
        {
            if(myid==0)
            {
                for(int k=1;k<mpi_num;k++)
                {
                   MPI_Send(velLocal_Collect[ip][kk].data(), cellNum[kk], MPI_DOUBLE, k, 666+kk, MPI_COMM_WORLD);
                }
            }
            else
            {
                MPI_Recv(velLocal_Collect[ip][kk].data(), cellNum[kk], MPI_DOUBLE, 0, 666+kk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}
MPI_Barrier(MPI_COMM_WORLD);

//npt is the number of IB points
for(label I=0; I<npt; I++)
{
    // Store the velocity components cal by interpolation
    scalar U_Interpolation_Tempx = 0.0;
    scalar U_Interpolation_Tempy = 0.0;

    //scalar P_Interpolation_Temp = 0.0; 
    //scalar P_Interpolation_Sum = 0.0;
    //std::vector<double> Temp3;
    //U_Interpolation.push_back(Temp3);

    // the number of cells near the IB points
    label Nsupport = meshSupport[I].size();
    //scalar Ntemp = meshSupport[I].size();

    for(label J=0; J<Nsupport; J++)
    {
        label N = meshSupport[I][J];
        label Np = meshSupportProcessor[I][J];

        U_Interpolation_Tempx += velLocal_Collect[0][Np][N]*Dirac[I][J]*h*h;
        U_Interpolation_Tempy += velLocal_Collect[1][Np][N]*Dirac[I][J]*h*h;
        //P_Interpolation_Temp += p[N]*Dirac[I][J]*h*h;
        //P_Interpolation_Sum += p[N];
    }
    //P_Interpolation_Temp = P_Interpolation_Sum/Nsupport;

    vector U_Interpolation_Temp(U_Interpolation_Tempx, U_Interpolation_Tempy, 0);
    U_Interpolation.push_back(U_Interpolation_Temp);

    //P_Interpolation.push_back(P_Interpolation_Temp);
}
/*
for(label I=0; I<nptZero; I++)
{
    label NN = meshSupport[0][I];
    label NP = meshSupportProcessor[0][I];
    Info << "No." << I+1 << ": " << NN << endl;
    Info << "Velocity: " << "Ux = " <<  velLocal_Collect[0][NP][NN] << "  " << "Uy = " << velLocal_Collect[1][NP][NN] << endl;
    //Info << "Pressure: " << p[NN] << endl;
}
*/
Info << "Lmarks[0] U_Interpolation: " << "U_Interpolation.x = " << U_Interpolation[0].x() << "  " << "U_Interpolation.y = " << U_Interpolation[0].y() << endl;
//Info << "Lmarks[0] P_Interpolation: " <<  P_Interpolation[0] << endl;

// store the IBpoints velocity and force(actually is acceleration)
std::vector<vector> F_Structure;
std::vector<vector> U_Structure;

Foam::dimensionedScalar dt = runTime.deltaT();
scalar dtt = dt.value();
//scalar dtt = 0.5;

Info << "Time step = " << dtt << endl;

for(label I=0; I<npt; I++)
{
    vector tempv(0.0, 0.0, 0.0);
    U_Structure.push_back(tempv);
}

Info << "Lmarks[0] U_Structure: " << "U_Structure.x = " << U_Structure[0].x() << "  " << "U_Structure.y = " << U_Structure[0].y() << endl;

for(label I=0; I<npt; I++)
{
    scalar F_Structurex = (U_Structure[I].x() - U_Interpolation[I].x())/dtt;
    scalar F_Structurey = (U_Structure[I].y() - U_Interpolation[I].y())/dtt;
    vector tempF(F_Structurex, F_Structurey, 0.0);
    F_Structure.push_back(tempF);
}

Info << "Lmarks[0] F_Structure: " << "F_Structurex = " << F_Structure[0].x() << "  " << "F_Structurey = " << F_Structure[0].y() << endl;
