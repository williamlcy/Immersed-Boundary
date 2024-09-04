// store the force of fluid
volVectorField F_Fluid
(
    IOobject
    (
        "F_Fluid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("F_Fluid", dimAcceleration, vector::zero),
    fixedValueFvPatchVectorField::typeName
);
//get the processor identifier
int myid=-1;
MPI_Comm_rank(MPI_COMM_WORLD,&myid);

//get the number of processors
int mpi_num=-1;
MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);

//get the mesh centers
const Foam::vectorField& meshCenters = mesh.cellCentres();

//get the number of cells in each processor
std::vector<int> cellNum(mpi_num);
//get the number of cells in the local processor
int local_cellnum=meshCenters.size();


//collect the number of cells in each processor
if(myid!=0)// if the present processor is not the root processor
{
    //send the number of cells in the present processor to the root processor
    MPI_Send(&local_cellnum, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
    /*
    local_cellnum: the adress of the variable to be sent
    1: the number of variables to be sent
    MPI_INT: the type of the variable to be sent (int)
    0: receive the variable in the root processor
    111: the tag of the message
    MPI_COMM_WORLD: the communicator
     */
}
else
{
    //collect the number of cells in the root processor
    cellNum[0]=local_cellnum;
    //collect the number of cells in the other processors
    for(int k=1;k<mpi_num;k++)
    {
     MPI_Recv(&cellNum[k], 1, MPI_INT, k, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

//send the number of cells in each processor to all the processors
if(myid==0)
{
    //root processor send
	for(int k=1;k<mpi_num;k++)
	{
	   MPI_Send(cellNum.data(), mpi_num, MPI_INT, k, 222, MPI_COMM_WORLD);
	}
}
else
{
    //other processors receive
	MPI_Recv(cellNum.data(), mpi_num, MPI_INT, 0, 222, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

//cal the total number of cells
int totolcellnum=0;
for(int k=0;k<mpi_num;k++)
{
   totolcellnum+=cellNum[k];
}


//collect the local processor mesh centers
std::vector<double> meshCenterLocal[3];
meshCenterLocal[0].resize(local_cellnum);
meshCenterLocal[1].resize(local_cellnum);
meshCenterLocal[2].resize(local_cellnum);
for(int j=0;j<local_cellnum;j++)
{
	meshCenterLocal[0][j]=meshCenters[j].x();
	meshCenterLocal[1][j]=meshCenters[j].y();
	meshCenterLocal[2][j]=meshCenters[j].z();
}

//collect the mesh centers of all the processors
std::vector<std::vector<double>> meshCenterLocal_Collect[3];

meshCenterLocal_Collect[0].resize(mpi_num);
meshCenterLocal_Collect[1].resize(mpi_num);
meshCenterLocal_Collect[2].resize(mpi_num);
for(int k=0;k<mpi_num;k++)
{
        meshCenterLocal_Collect[0][k].resize(cellNum[k]);
	meshCenterLocal_Collect[1][k].resize(cellNum[k]);
        meshCenterLocal_Collect[2][k].resize(cellNum[k]);
}

//make sure all the processors have finished the previous operations
//MPI_COMM_WORLD: the communicator(include all the processors
MPI_Barrier(MPI_COMM_WORLD);

//operation and exchange mesh center coordination for three directions
//ip: the direction of the mesh center(x,y,z)
for(int ip=0;ip<3;ip++)
{
	if(myid!=0)
	{
		MPI_Send(meshCenterLocal[ip].data(), local_cellnum, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);
	}
	else
	{
		meshCenterLocal_Collect[ip][0]=meshCenterLocal[ip];
		for(int k=1;k<mpi_num;k++)
		{
		 MPI_Recv(meshCenterLocal_Collect[ip][k].data(), cellNum[k], MPI_DOUBLE, k, 333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(int kk=0;kk<mpi_num;kk++)
	{
		if(myid==0)
		{
			for(int k=1;k<mpi_num;k++)
			{
			   MPI_Send(meshCenterLocal_Collect[ip][kk].data(), cellNum[kk], MPI_DOUBLE, k, 444+kk, MPI_COMM_WORLD);
			}
		}
		else
		{
			MPI_Recv(meshCenterLocal_Collect[ip][kk].data(), cellNum[kk], MPI_DOUBLE, 0, 444+kk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);

