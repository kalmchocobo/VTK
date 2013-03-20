/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSignedImplicitModeller.cxx

  Copyright (c) Alex Mansfield
  All rights reserved.
  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
#include "vtkSignedImplicitModeller.h"

#include "vtkBox.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkClipPolyData.h"
#include "vtkCommand.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkImageData.h"
#include "vtkImageIterator.h"
#include "vtkImageProgressIterator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMath.h"
#include "vtkMultiThreader.h"
#include "vtkMutexLock.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPolyDataNormals.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"

#include <math.h>
#include <limits>

vtkStandardNewMacro(vtkSignedImplicitModeller);

struct vtkSignedImplicitModellerAppendInfo
{
	vtkSignedImplicitModeller  *Modeller;
	vtkDataSet                **Input;
	double                      MaximumDistance;
};

//----------------------------------------------------------------------------
// Construct with sample dimensions=(50,50,50), and so that model bounds are
// automatically computed from the input. Capping is turned on with CapValue
// equal to a large positive number.
vtkSignedImplicitModeller::vtkSignedImplicitModeller()
{
	this->MaximumDistance = 0.1;

	this->ModelBounds[0] = 0.0;
	this->ModelBounds[1] = 0.0;
	this->ModelBounds[2] = 0.0;
	this->ModelBounds[3] = 0.0;
	this->ModelBounds[4] = 0.0;
	this->ModelBounds[5] = 0.0;
	this->BoundsComputed = 0;

	this->Spacing = 0;
	this->SampleDimensions[0] = 50;
	this->SampleDimensions[1] = 50;
	this->SampleDimensions[2] = 50;

	this->Capping = 1;
	this->OutputScalarType = VTK_FLOAT;
	this->CapValue = this->GetScalarTypeMax( this->OutputScalarType );
	this->ScaleToMaximumDistance = 0; // only used for non-float output type
  
	this->DataAppended = 0;
	this->AdjustBounds = 1;
	this->AdjustDistance = 0.0125;
	this->AdjustBoundsAbsolute = 0;
	this->AdjustDistanceAbsolute = 0.0;
	this->BoxesSegmentation.clear();

	// this->ProcessMode = VTK_CELL_MODE;
	this->LocatorMaxLevel = 5;

	this->Threader        = vtkMultiThreader::New();
	this->NumberOfThreads = this->Threader->GetNumberOfThreads();

}

vtkSignedImplicitModeller::~vtkSignedImplicitModeller()
{
	if (this->Threader)
    {
		this->Threader->Delete();
    }
}

int vtkSignedImplicitModeller::ComputeSignedDistance(vtkPolyData * input, int cellNum, Barycentric const & barycentric, double x[3], double & sdist)
{
	// initialize outputs
	sdist = std::numeric_limits<double>::max();

	// first find the closest point
	double cp[3], cp_barycentric[3];

	// get cell
	vtkCell * cell = input->GetCell(cellNum);

	// get plane normal
	double planeNormal[3];
	input->GetCellData()->GetNormals()->GetTuple(cellNum,planeNormal);
	
	// get the 3 vertices
	double p0[3], p1[3], p2[3];
	cell->GetPoints()->GetPoint(0, p0);
	cell->GetPoints()->GetPoint(1, p1);
	cell->GetPoints()->GetPoint(2, p2);

	/* check for case that closest point is inside the triangle */

	// project point x to plane, closest point is cp
	double x_p[3];
	vtkPlane::GeneralizedProjectPoint(x, p2, planeNormal, x_p);

	// compute Barycentric co-ordinates
	double x_p_b[3];
	barycentric.compute(x_p,x_p_b);
	
	if ( x_p_b[0] >= 0.0 and x_p_b[0] <= 1.0 and
		 x_p_b[1] >= 0.0 and x_p_b[1] <= 1.0 and
		 x_p_b[2] >= 0.0 and x_p_b[2] <= 1.0 )
    {
		// closest point in the triangle is x_plane
		for (int i=0; i<3; i++)
		{
			cp[i]             = x_p[i];
			cp_barycentric[i] = x_p_b[i];
		}
	}
	else
	{
		/* otherwise, closest point lies on boundary of triangle */

		double dist[3], t[3], cp_l[9];

		dist[0] = vtkLine::DistanceToLine(x,p2,p1,t[0],&(cp_l[0]));
		dist[1] = vtkLine::DistanceToLine(x,p2,p0,t[1],&(cp_l[3]));
		dist[2] = vtkLine::DistanceToLine(x,p1,p0,t[2],&(cp_l[6]));

		// find closest line
		int closest = 0;
		if(dist[1] < dist[0])
			closest = 1;
		if(dist[2] < dist[1] and dist[2] < dist[0])
			closest = 2;

		for (int i=0; i<3; i++)
			cp[i] = cp_l[3*closest + i];

		barycentric.compute(cp, cp_barycentric);
	}

	// we now have cp and cp_barycentric

	// compute distance from x to cp
	sdist = std::sqrt(vtkMath::Distance2BetweenPoints(x,cp));

	// compute corresponding normal
	double normal[3];

	// get the 3 vertex normals
	// double vNormal[9];
	// for(int i=0; i<3; i++)
	// 	input->GetPointData()->GetNormals()->GetTuple(cell->GetPointIds()->GetId(i),&(vNormal[3*i]));

	// for each dimension
	for(int i=0; i<3; i++)
	{
		normal[i] = planeNormal[i];

		// normal[i] = 0.0;

		// // for each vertex of the triangle
		// for(int j=0; j<3; j++)
		// 	normal[i] += cp_barycentric[j] * vNormal[3*j + i];
	}

	// compute sign
	double diffvec[3];
	vtkMath::Subtract(x,cp,diffvec);

	double dot_product = vtkMath::Dot(diffvec,normal);
	
	// negative distance if dot_product<0
	// but only if can be confident i.e. abs(dot_product) not too small
	if( (dot_product<0.0) and (1e-12 < (dot_product)*(dot_product)) )
		sdist *= -1.0;
	
	//std::printf("proposed_distance %f dot_product %f bcoords %f %f %f vec %f %f %f normal %f %f %f cp %f %f %f \n",sdist,vtkMath::Dot(diffvec,normal),cp_barycentric[0],cp_barycentric[1],cp_barycentric[2],diffvec[0],diffvec[1],diffvec[2],normal[0],normal[1],normal[2],cp[0],cp[1],cp[2]);
	
	return 1;
}

void vtkSignedImplicitModeller::SetOutputScalarType(int type)
{ 
	double scalarMax;

	vtkDebugMacro(<< this->GetClassName() << " (" << this <<
				  "): setting OutputScalarType to " << type); 

	scalarMax = this->GetScalarTypeMax(type);
	if (scalarMax) // legal type
    {
		int modified = 0;
		if (this->CapValue != scalarMax)
		{
			this->CapValue = scalarMax;
			modified = 1;
		}
		if (this->OutputScalarType != type) 
		{ 
			this->OutputScalarType = type; 
			modified = 1;
		}
		if (modified)
		{
			this->Modified(); 
		}
    } 
} 

void vtkSignedImplicitModeller::SetCapValue(double value)
{
	vtkDebugMacro(<< this->GetClassName() << " (" << this <<
				  "): setting CapValue to " << value);
	// clamp to between 0 and max for scalar type
	double max = this->GetScalarTypeMax(this->OutputScalarType);
	if (this->CapValue != (value < 0 ? 0 : (value > max ? max : value)))
    {
		this->CapValue = (value < 0 ? 0 : (value > max ? max : value));
		this->Modified();
    }
}


double vtkSignedImplicitModeller::GetScalarTypeMax(int type)
{
	switch (type)
    {
    case VTK_UNSIGNED_CHAR:  return (double)VTK_UNSIGNED_CHAR_MAX;
    case VTK_CHAR:           return (double)VTK_CHAR_MAX;
    case VTK_UNSIGNED_SHORT: return (double)VTK_UNSIGNED_SHORT_MAX;
    case VTK_SHORT:          return (double)VTK_SHORT_MAX;
    case VTK_UNSIGNED_INT:   return (double)VTK_UNSIGNED_INT_MAX;
    case VTK_INT:            return (double)VTK_INT_MAX;
    case VTK_UNSIGNED_LONG:  return (double)VTK_UNSIGNED_LONG_MAX;
    case VTK_LONG:           return (double)VTK_LONG_MAX;
    case VTK_FLOAT:          return (double)VTK_FLOAT_MAX;
    case VTK_DOUBLE:         return (double)VTK_DOUBLE_MAX;
    default: return 0;
    }
}

void vtkSignedImplicitModeller::SetBoxesSegmentation(std::list<std::vector<double> > const & _BoxesSegmentation)
{
	BoxesSegmentation = _BoxesSegmentation;
	return;
}

void vtkSignedImplicitModeller::GetBoxesSegmentation(std::list<std::vector<double> > & _BoxesSegmentation)
{
	_BoxesSegmentation = BoxesSegmentation;
	return;
}

//----------------------------------------------------------------------------
// Initialize the filter for appending data. You must invoke the
// StartAppend() method before doing successive Appends(). It's also a
// good idea to manually specify the model bounds; otherwise the input
// bounds for the data will be used.
void vtkSignedImplicitModeller::StartAppend()
{
	this->StartAppend(0);
}

void vtkSignedImplicitModeller::StartAppend(int internal)
{
	vtkIdType numPts;
	vtkIdType i;
	double maxDistance;

	if (!internal)
    {
		// we must call update information because we can't be sure that 
		// it has been called.
		this->UpdateInformation();
    }
	vtkStreamingDemandDrivenPipeline::SetUpdateExtent(this->GetOutputInformation(0), vtkStreamingDemandDrivenPipeline::GetWholeExtent(this->GetOutputInformation(0)));
  
	vtkDebugMacro(<< "Initializing data");
	this->AllocateOutputData(this->GetOutput(), this->GetOutputInformation(0));
	this->UpdateProgress(0.0);
	this->DataAppended = 1;

	numPts = this->SampleDimensions[0] * this->SampleDimensions[1] * this->SampleDimensions[2];

	// initialize output to CapValue at each location
	maxDistance = this->CapValue;
	vtkDataArray *newScalars = this->GetOutput()->GetPointData()->GetScalars();
	for (i=0; i<numPts; i++)
    {
		newScalars->SetComponent(i, 0, maxDistance);
    }
	
}

template <class OT>
void vtkSignedImplicitModeller::SetOutputDistance(double distance, OT *outputValue, double capValue, double scaleFactor)
{
	// for now, just doing "normal" cast... could consider doing round?
	if (scaleFactor)  // need to scale the distance
    {
		*outputValue = static_cast<OT>(distance * scaleFactor);
    }
	else 
    {
		if (capValue and std::fabs(distance) > capValue) // clamping iff non-float type
			distance = (distance>0 - distance<0) * capValue;

		*outputValue = static_cast<OT>(distance);
    }
}

// Convert distance as stored in output (could be scaled and/or non-double type) to double distance with correct scaling
template <class OT>
void vtkSignedImplicitModeller::GetOutputDistance(OT outputValue, double & distance, double & distance2, double scaleFactor)
{
	distance = (double)outputValue;

	if (scaleFactor) 
		distance *= scaleFactor;

	distance2 = distance * distance;
}


//----------------------------------------------------------------------------
// Templated append for VTK_CELL_MODE process mode and any type of output data
template <class OT>
void vtkSignedImplicitModellerAppendExecute(vtkSignedImplicitModeller *self,
											vtkDataSet *input, vtkImageData *outData,
											double maxDistance, OT *)
{
	double const * const spacing = outData->GetSpacing();
	double const * const origin = outData->GetOrigin();    
	int    const * const sampleDimensions = self->GetSampleDimensions();

	double const maxDistance2 = maxDistance * maxDistance;    
	
	// so we know how to scale if desired
	double scaleFactor = 0; // 0 used to indicate not scaling
	double toDoubleScaleFactor = 0; // 0 used to indicate not scaling 
	double capValue = 0; // 0 used to indicate not clamping (float or double)
	if (self->GetOutputScalarType() != VTK_FLOAT and self->GetOutputScalarType() != VTK_DOUBLE)
    {
		capValue = self->GetCapValue();
		if (self->GetScaleToMaximumDistance())
		{
			scaleFactor = capValue / maxDistance;
			toDoubleScaleFactor = maxDistance / capValue;
		}
    }

	vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(input);

	// get BoxesSegmentation
	std::list<std::vector<double> > BoxesSegmentation;
	self->GetBoxesSegmentation(BoxesSegmentation);
	bool const UseSegmentation = BoxesSegmentation.size()>0;

	// remove segmentation boxes
	if( UseSegmentation )
	{
		std::list<std::vector<double> >::iterator boxes_it = BoxesSegmentation.begin();
		for(; boxes_it != BoxesSegmentation.end(); boxes_it++)
		{
			// create bounding box rectangle
			vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
			box->SetBounds(boxes_it->data());

			vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
			clipper->SetInputData(poly);
			clipper->SetClipFunction(box);
			clipper->Update();
			poly = vtkPolyData::SafeDownCast(clipper->GetOutput());
		}
	}

	// compute face and vertex normals if we don't have them already
	vtkSmartPointer<vtkPolyDataNormals> normalsGen = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalsGen->SetInputData(poly);
	normalsGen->ComputeCellNormalsOn();
	normalsGen->ComputePointNormalsOff();

	// I seem to get the best results with all these fancy options turned off
	normalsGen->SplittingOff();
	normalsGen->ConsistencyOff();
	normalsGen->AutoOrientNormalsOff();
	normalsGen->NonManifoldTraversalOff();

	// generate normals
	normalsGen->Update();
	poly = vtkPolyData::SafeDownCast(normalsGen->GetOutput());

	// obtain point/face mappings ***
	// poly->BuildLinks();

	//
	// Traverse all cells; computing distance function on volume points.
	//

	int const updateTime = std::max(1, (int)((double)poly->GetNumberOfCells() / (double)50.0));  // update every 2%

	// for every cell (triangle in the mesh)
	for (vtkIdType cellNum = 0; cellNum < poly->GetNumberOfCells(); cellNum++)
    {
		// compute volume of space in which this triangle could be a closest point
		// (maxDistance is the maximum distance in which SDF calculated)
		vtkCell * const cell = poly->GetCell(cellNum);
	
		// get the 3 vertices
		double p0[3], p1[3], p2[3];
		cell->GetPoints()->GetPoint(0, p0);
		cell->GetPoints()->GetPoint(1, p1);
		cell->GetPoints()->GetPoint(2, p2);

		// create Barycentric co-ordinates
		vtkSignedImplicitModeller::Barycentric barycentric(p0, p1, p2);

		if(barycentric.get_initialized())
		{
			double adjBounds[6];
			double * const bounds = cell->GetBounds();
			for (int i=0; i<3; i++)
			{
				adjBounds[2*i] = bounds[2*i] - maxDistance;
				adjBounds[2*i+1] = bounds[2*i+1] + maxDistance;
			}

			// compute range of voxels to test
			int outExt[6];
			for (int i=0; i<3; i++)
			{
				outExt[i*2]     = (int) ((double)(adjBounds[2*i]   - origin[i]) / spacing[i]);
				outExt[i*2 + 1] = (int) ((double)(adjBounds[2*i+1] - origin[i]) / spacing[i]);
				if (outExt[i*2] < 0)
				{
					outExt[i*2] = 0;
				}
				if (outExt[i*2 + 1] >= sampleDimensions[i])
				{
					outExt[i*2 + 1] = sampleDimensions[i] - 1;
				}
			}

			// create iterator over SDF image outData in range of voxels to test outExt
			vtkImageIterator<OT> outIt(outData, outExt);

			// loop through all voxels to test
			// compute x, the centre of the voxel
			double x[3];
			for (int k = outExt[4]; k <= outExt[5]; k++) 
			{
				x[2] = spacing[2] * k + origin[2];
				for (int j = outExt[2]; j <= outExt[3]; j++)
				{
					x[1] = spacing[1] * j + origin[1];
					OT* outSI = outIt.BeginSpan();
					for (int i = outExt[0]; i <= outExt[1]; i++) 
					{
						x[0] = spacing[0] * i + origin[0];

						// if(i==20 and j==14 and k==8) //***
						// {
						// get current distance and distance2 at this voxel
						double prevDistance, prevDistance2;
						self->GetOutputDistance(*outSI, prevDistance, prevDistance2, toDoubleScaleFactor);

						//std::printf("(%i %i %i) (%f %f %f) min_distance %f ",i,j,k,x[0],x[1],x[2],distance);

						// compute distance to closest point on the cell from x
						// if computation occurs without numerical error AND if distance less than maximum distance
						double distance;
						if(self->ComputeSignedDistance(poly, cellNum, barycentric, x, distance)>0 and
							distance*distance <= maxDistance2)
						{
							// deal with thin structures
							// if distance same but signs different, assume positive sign (outside)
							bool thin_structure = false;
							if( (((distance*distance)-prevDistance2)*((distance*distance)-prevDistance2)) < 1e-12 and
								( (distance>0 and prevDistance<0) or (distance<0 and prevDistance>0) ) )
							{
								thin_structure = true;

								// choose sign +ve
								if(distance<0)
									distance *= -1.0;
							}

							// set distance at this voxel in SDF image
							if(distance*distance < prevDistance2 or thin_structure)
								self->SetOutputDistance(distance, outSI, capValue, scaleFactor);
						}
						
						// }
						outSI++;
					}
					outIt.NextSpan();
				}
			}
		}

		if (cellNum % updateTime == 0)
			self->UpdateProgress(double(cellNum + 1) / poly->GetNumberOfCells());
	}
}

// Append a data set to the existing output. To use this function,
// you'll have to invoke the StartAppend() method before doing
// successive appends. It's also a good idea to specify the model
// bounds; otherwise the input model bounds is used. When you've
// finished appending, use the EndAppend() method.
void vtkSignedImplicitModeller::Append(vtkDataSet *input)
{
	vtkDebugMacro(<< "Appending data");

	vtkImageData *output = this->GetOutput();

	if ( !this->BoundsComputed )
	{
		this->ComputeModelBounds(input);
	}

	// if (this->ProcessMode == VTK_CELL_MODE)
	// {
	if (!output->GetPointData()->GetScalars())
	{
		vtkErrorMacro("Sanity check failed.");
		return;
	}
 
	switch (this->OutputScalarType)
	{
		vtkTemplateMacro(
			vtkSignedImplicitModellerAppendExecute( this, 
													input, 
													output, 
													this->InternalMaxDistance, 
													static_cast<VTK_TT *>(0)));
	}
    // }
	// else
    // {
	// 	vtkSignedImplicitModellerAppendInfo info;
	// 	double minZ, maxZ;
	// 	int slabMin, slabMax, slabSize, i;
	// 	vtkClipPolyData **minClipper = NULL, **maxClipper = NULL; 
	// 	vtkPlane ** minPlane = NULL, **maxPlane = NULL;
	// 	double *spacing, *origin;

	// 	spacing = output->GetSpacing();
	// 	origin = output->GetOrigin();
    
	// 	// Use a MultiThreader here, splitting the volume into slabs to be processed
	// 	// by the separate threads
    
	// 	// Set the number of threads to use,
	// 	// then set the execution method and do it.
	// 	this->Threader->SetNumberOfThreads( this->NumberOfThreads );
    
	// 	// set up the info object for the thread
	// 	info.Modeller = this;
	// 	info.MaximumDistance = this->InternalMaxDistance;

	// 	info.Input = new vtkDataSet* [this->NumberOfThreads];
	// 	if (this->NumberOfThreads == 1)
	// 	{
	// 		info.Input[0] = input;
	// 	}
	// 	else
	// 	{
	// 		// if not PolyData, then copy the input for each thread
	// 		if ( input->GetDataObjectType() != VTK_POLY_DATA )
	// 		{
	// 			for (i = 0; i < this->NumberOfThreads; i++)
	// 			{
	// 				switch( input->GetDataObjectType() )
	// 				{
	// 				case VTK_STRUCTURED_GRID:
	// 					info.Input[i] = vtkStructuredGrid::New();
	// 					break;
	// 				case VTK_IMAGE_DATA:
	// 					info.Input[i] = vtkImageData::New();
	// 					break;
	// 				case VTK_UNSTRUCTURED_GRID:
	// 					info.Input[i] = vtkUnstructuredGrid::New();
	// 					break;
	// 				case VTK_RECTILINEAR_GRID:
	// 					info.Input[i] = vtkRectilinearGrid::New();
	// 					break;
	// 				default:
	// 					vtkErrorMacro(<<"Unexpected DataSet type!");
	// 					return;
	// 				}
	// 				info.Input[i]->CopyStructure(input);
	// 			}
	// 		}
	// 		else // break up the input data into slabs to help ensure thread safety

	// 		{
	// 			minClipper = new vtkClipPolyData* [this->NumberOfThreads];
	// 			maxClipper = new vtkClipPolyData* [this->NumberOfThreads];
	// 			minPlane = new vtkPlane* [this->NumberOfThreads];
	// 			maxPlane = new vtkPlane* [this->NumberOfThreads];
        
	// 			slabSize = this->SampleDimensions[2] / this->NumberOfThreads;
	// 			if (slabSize == 0) // in case threadCount >  SampleDimensions[2]
	// 			{
	// 				slabSize = 1;
	// 			}

	// 			for (i = 0; i < this->NumberOfThreads; i++)
	// 			{
	// 				minPlane[i] = maxPlane[i] = 0;
	// 				//////////////////////////////////////////////////
	// 				// do the 1st clip
	// 				slabMin = i * slabSize;
	// 				if (slabMin >= this->SampleDimensions[2])
	// 				{
	// 					break;
	// 				}

	// 				// get/clip input cells in this slab + maxDistance+ 
	// 				minZ = spacing[2] * slabMin + origin[2] - this->InternalMaxDistance*1.00001;
	// 				if (minZ < this->ModelBounds[4])
	// 				{
	// 					minZ = this->ModelBounds[4];
	// 				}

	// 				minPlane[i] = vtkPlane::New();
	// 				minPlane[i]->SetNormal(0.0, 0.0, -1.0);
	// 				minPlane[i]->SetOrigin(0.0, 0.0, minZ);

	// 				minClipper[i] = vtkClipPolyData::New();
	// 				minClipper[i]->SetInputData((vtkPolyData *)input);
	// 				minClipper[i]->SetClipFunction(minPlane[i]);
	// 				minClipper[i]->SetValue( 0.0 );
	// 				minClipper[i]->InsideOutOn();
	// 				minClipper[i]->Update();

	// 				if ( minClipper[i]->GetOutput()->GetNumberOfCells() == 0 )
	// 				{
	// 					info.Input[i] = NULL;
	// 					maxPlane[i] = NULL;
	// 					continue;
	// 				}
	// 				minClipper[i]->ReleaseDataFlagOn();

	// 				//////////////////////////////////////////////////
	// 				// do the 2nd clip
	// 				slabMax = slabMin + slabSize - 1;
	// 				if (i == this->NumberOfThreads - 1)
	// 				{
	// 					slabMax = this->SampleDimensions[2] - 1;
	// 				}
          
	// 				maxZ = spacing[2] * slabMax + origin[2] + this->InternalMaxDistance*1.00001;
	// 				if (maxZ > this->ModelBounds[5])
	// 				{
	// 					maxZ = this->ModelBounds[5];
	// 				}
	// 				maxPlane[i] = vtkPlane::New();
	// 				maxPlane[i]->SetNormal(0.0, 0.0, 1.0);
	// 				maxPlane[i]->SetOrigin(0.0, 0.0, maxZ);

	// 				maxClipper[i] = vtkClipPolyData::New();
	// 				maxClipper[i]->SetInputData(minClipper[i]->GetOutput());
	// 				maxClipper[i]->SetClipFunction(maxPlane[i]);
	// 				maxClipper[i]->SetValue( 0.0 );
	// 				maxClipper[i]->InsideOutOn();
	// 				maxClipper[i]->Update();

	// 				if ( maxClipper[i]->GetOutput()->GetNumberOfCells() == 0 )
	// 				{
	// 					info.Input[i] = NULL;
	// 				}
	// 				else
	// 				{
	// 					info.Input[i] = maxClipper[i]->GetOutput();
	// 				}
	// 			}
	// 			for (i++; i < this->NumberOfThreads; i++)
	// 			{
	// 				minPlane[i] = maxPlane[i] = 0;
	// 			}
	// 		}
	// 	}

	// 	this->Threader->SetSingleMethod( vtkSignedImplicitModeller_ThreadedAppend, 
	// 									 (void *)&info);
	// 	this->Threader->SingleMethodExecute();

	// 	// cleanup
	// 	if (this->NumberOfThreads > 1)
	// 	{
	// 		if ( input->GetDataObjectType() != VTK_POLY_DATA )
	// 		{
	// 			for (i = 0; i < this->NumberOfThreads; i++)
	// 			{
	// 				info.Input[i]->Delete();
	// 			}
	// 		}
	// 		else
	// 		{
	// 			for (i = 0; i < this->NumberOfThreads; i++)
	// 			{
	// 				if (minPlane[i])
	// 				{
	// 					minPlane[i]->Delete();
	// 					minClipper[i]->Delete();
	// 				}
	// 				if (maxPlane[i])
	// 				{
	// 					maxPlane[i]->Delete();
	// 					maxClipper[i]->Delete();
	// 				}
	// 			}
	// 			delete [] minPlane;
	// 			delete [] maxPlane;
	// 			delete [] minClipper;
	// 			delete [] maxClipper;
	// 		}
	// 	}
	// 	delete [] info.Input;
    // }
}

//----------------------------------------------------------------------------
// Method completes the append process (does the capping if requested).
void vtkSignedImplicitModeller::EndAppend()
{
	vtkDataArray *newScalars;
	vtkDebugMacro(<< "End append");
  
	if (!(newScalars =this->GetOutput()->GetPointData()->GetScalars()))
    {
		vtkErrorMacro("Sanity check failed.");
		return;
    }

	if ( this->Capping )
    {
		this->Cap(newScalars);
    }
	this->UpdateProgress(1.0);
}

//----------------------------------------------------------------------------
int vtkSignedImplicitModeller::RequestInformation (
	vtkInformation * vtkNotUsed(request),
	vtkInformationVector ** vtkNotUsed( inputVector ),
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	vtkDataObject::SetPointDataActiveScalarInfo(outInfo, this->OutputScalarType, 1);

	// Set volume origin
	double origin[3];
	for(int i=0; i<3; i++)
		origin[i] = this->ModelBounds[2*i];

	outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);

	// Set data spacing
	double vSpacing[3];
	if(this->Spacing > 0)
	{
		for (int i=0; i<3; i++)
			vSpacing[i] = this->Spacing;
	} else {
		for (int i=0; i<3; i++)
			vSpacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i]) / (this->SampleDimensions[i] - 1);
	}

	outInfo->Set(vtkDataObject::SPACING(),vSpacing,3);

	// Set extent
	if(this->Spacing > 0)
	{
		for (int i=0; i<3; i++)
			this->SampleDimensions[i] = std::ceil(0.5+((this->ModelBounds[2*i+1]-this->ModelBounds[2*i])/vSpacing[i]));
	}

	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
				 0, this->SampleDimensions[0]-1,
				 0, this->SampleDimensions[1]-1,
				 0, this->SampleDimensions[2]-1);

	return 1;
}

//----------------------------------------------------------------------------
int vtkSignedImplicitModeller::RequestData(
	vtkInformation* vtkNotUsed( request ),
	vtkInformationVector** inputVector,
	vtkInformationVector* vtkNotUsed( outputVector ))
{
	// get the input
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
  
	vtkDebugMacro(<< "Executing implicit model");

	if (input == NULL)
    {
		// we do not want to release the data because user might
		// have called Append ...
		return 0;
    }

	this->StartAppend(1);
	this->Append(input);
	this->EndAppend();

	return 1;
}

// Compute ModelBounds from input geometry.
double vtkSignedImplicitModeller::ComputeModelBounds(vtkDataSet *input)
{
	double *bounds, maxDist;
	int i;
	vtkImageData *output = this->GetOutput();
  
	// compute model bounds if not set previously
	if ( this->ModelBounds[0] >= this->ModelBounds[1] ||
		 this->ModelBounds[2] >= this->ModelBounds[3] ||
		 this->ModelBounds[4] >= this->ModelBounds[5] )
    {
		if (input != NULL)
		{
			bounds = input->GetBounds();
		}
		else
		{
			vtkDataSet *dsInput = vtkDataSet::SafeDownCast(this->GetInput());
			if (dsInput != NULL)
			{
				bounds = dsInput->GetBounds();
			}
			else
			{
				vtkErrorMacro( 
					<< "An input must be specified to Compute the model bounds.");
				return VTK_FLOAT_MAX;
			}
		}
    }
	else
    {
		bounds = this->ModelBounds;
    }

	for (maxDist=0.0, i=0; i<3; i++)
    {
		if ( (bounds[2*i+1] - bounds[2*i]) > maxDist )
		{
			maxDist = bounds[2*i+1] - bounds[2*i];
		}
    }

	for (i=0; i<3; i++)
	{
		this->ModelBounds[2*i] = bounds[2*i];
		this->ModelBounds[2*i+1] = bounds[2*i+1];
	}

	// adjust bounds so model fits strictly inside (only if not set previously)
	if ( this->AdjustBounds )
    {
		for (i=0; i<3; i++)
		{
			this->ModelBounds[2*i] -= maxDist*this->AdjustDistance;
			this->ModelBounds[2*i+1] += maxDist*this->AdjustDistance;
		}
    }

	if ( this->AdjustBoundsAbsolute )
    {
		for (i=0; i<3; i++)
		{
			this->ModelBounds[2*i] -= this->AdjustDistanceAbsolute;
			this->ModelBounds[2*i+1] += this->AdjustDistanceAbsolute;
		}
    }

	maxDist *= this->MaximumDistance;

	vtkInformation *outInfo = this->GetExecutive()->GetOutputInformation(0);

	// Set volume origin
	output->SetOrigin(this->ModelBounds[0], this->ModelBounds[2], this->ModelBounds[4]);
	outInfo->Set(vtkDataObject::ORIGIN(), this->ModelBounds[0], this->ModelBounds[2], this->ModelBounds[4]);

	// Set data spacing
	double vSpacing[3];
	if(this->Spacing > 0)
	{
		for (i=0; i<3; i++)
			vSpacing[i] = this->Spacing;
	} else {
		for (i=0; i<3; i++)
			vSpacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i]) / (this->SampleDimensions[i] - 1);
	}

	output->SetSpacing(vSpacing);
	outInfo->Set(vtkDataObject::SPACING(),vSpacing,3);

	// Compute SampleDimensions
	if(this->Spacing > 0)
	{
		for (i=0; i<3; i++)
			this->SampleDimensions[i] = std::ceil(0.5+((this->ModelBounds[2*i+1]-this->ModelBounds[2*i])/vSpacing[i]));
	}

	output->SetExtent(0,this->SampleDimensions[0]-1,0,this->SampleDimensions[1]-1,0,this->SampleDimensions[2]-1);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
				 0, this->SampleDimensions[0]-1,
				 0, this->SampleDimensions[1]-1,
				 0, this->SampleDimensions[2]-1);

	this->BoundsComputed = 1;
	this->InternalMaxDistance = maxDist;

	return maxDist;  
}

//----------------------------------------------------------------------------
// Set the i-j-k dimensions on which to sample the distance function.
void vtkSignedImplicitModeller::SetSampleDimensions(int i, int j, int k)
{
	int dim[3];

	dim[0] = i;
	dim[1] = j;
	dim[2] = k;

	this->SetSampleDimensions(dim);
}

//----------------------------------------------------------------------------
void vtkSignedImplicitModeller::SetSampleDimensions(int dim[3])
{
	int dataDim, i;

	vtkDebugMacro(<< " setting SampleDimensions to (" << dim[0] << "," << dim[1] << "," << dim[2] << ")");

	if ( dim[0] != this->SampleDimensions[0] ||
		 dim[1] != this->SampleDimensions[1] ||
		 dim[2] != this->SampleDimensions[2] )
    {
		if ( dim[0]<1 || dim[1]<1 || dim[2]<1 )
		{
			vtkErrorMacro (<< "Bad Sample Dimensions, retaining previous values");
			return;
		}

		for (dataDim=0, i=0; i<3 ; i++)
		{
			if (dim[i] > 1)
			{
				dataDim++;
			}
		}

		if ( dataDim  < 3 )
		{
			vtkErrorMacro(<<"Sample dimensions must define a volume!");
			return;
		}

		for ( i=0; i<3; i++)
		{
			this->SampleDimensions[i] = dim[i];
		}

		this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkSignedImplicitModeller::Cap(vtkDataArray *s)
{
	int i,j,k;
	int idx;
	int d01=this->SampleDimensions[0]*this->SampleDimensions[1];

// i-j planes
	k = 0;
	for (j=0; j<this->SampleDimensions[1]; j++)
    {
		for (i=0; i<this->SampleDimensions[0]; i++)
		{
			s->SetComponent(i+j*this->SampleDimensions[0],0, this->CapValue);
		}
    }
	k = this->SampleDimensions[2] - 1;
	idx = k*d01;
	for (j=0; j<this->SampleDimensions[1]; j++)
    {
		for (i=0; i<this->SampleDimensions[0]; i++)
		{
			s->SetComponent(idx+i+j*this->SampleDimensions[0], 0, this->CapValue);
		}
    }
	// j-k planes
	i = 0;
	for (k=0; k<this->SampleDimensions[2]; k++)
    {
		for (j=0; j<this->SampleDimensions[1]; j++)
		{
			s->SetComponent(j*this->SampleDimensions[0]+k*d01,0,this->CapValue);
		}
    }
	i = this->SampleDimensions[0] - 1;
	for (k=0; k<this->SampleDimensions[2]; k++)
    {
		for (j=0; j<this->SampleDimensions[1]; j++)
		{
			s->SetComponent(i+j*this->SampleDimensions[0]+k*d01,0, this->CapValue);
		}
    }
	// i-k planes
	j = 0;
	for (k=0; k<this->SampleDimensions[2]; k++)
    {
		for (i=0; i<this->SampleDimensions[0]; i++)
		{
			s->SetComponent(i+k*d01,0, this->CapValue);
		}
    }
	j = this->SampleDimensions[1] - 1;
	idx = j*this->SampleDimensions[0];
	for (k=0; k<this->SampleDimensions[2]; k++)
    {
		for (i=0; i<this->SampleDimensions[0]; i++)
		{
			s->SetComponent(idx+i+k*d01,0, this->CapValue);
		}
    }
}

//----------------------------------------------------------------------------
int vtkSignedImplicitModeller::FillInputPortInformation(
	int vtkNotUsed( port ), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
	info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
	return 1;
}

//----------------------------------------------------------------------------
int vtkSignedImplicitModeller::ProcessRequest(vtkInformation* request,
											  vtkInformationVector** inputVector,
											  vtkInformationVector* outputVector)
{
	// If we have no input then we will not generate the output because
	// the user already called StartAppend/Append/EndAppend.
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED()))
    {
		if(inputVector[0]->GetNumberOfInformationObjects() == 0)
		{
			vtkInformation* outInfo = outputVector->GetInformationObject(0);
			outInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
		}
		return 1;
    }
	else if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
		if(inputVector[0]->GetNumberOfInformationObjects() == 0)
		{
			return 1;
		}
    }
	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

void vtkSignedImplicitModeller::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	os << indent << "Maximum Distance: " << this->MaximumDistance << "\n";
	os << indent << "OutputScalarType: " << this->OutputScalarType << "\n";
	os << indent << "Sample Dimensions: (" << this->SampleDimensions[0] << ", "
	   << this->SampleDimensions[1] << ", "
	   << this->SampleDimensions[2] << ")\n";
	os << indent << "ModelBounds: \n";
	os << indent << "  Xmin,Xmax: (" << this->ModelBounds[0] << ", " 
	   << this->ModelBounds[1] << ")\n";
	os << indent << "  Ymin,Ymax: (" << this->ModelBounds[2] << ", " 
	   << this->ModelBounds[3] << ")\n";
	os << indent << "  Zmin,Zmax: (" << this->ModelBounds[4] << ", " 
	   << this->ModelBounds[5] << ")\n";

	os << indent << "ScaleToMaximumDistance: " << (this->ScaleToMaximumDistance ? "On\n" : "Off\n");
	os << indent << "AdjustBounds: " << (this->AdjustBounds ? "On\n" : "Off\n");
	os << indent << "AdjustDistance: " << this->AdjustDistance << "\n";
	os << indent << "AdjustBoundsAbsolute: " << (this->AdjustBoundsAbsolute ? "On\n" : "Off\n");
	os << indent << "AdjustDistanceAbsolute: " << this->AdjustDistanceAbsolute << "\n";
	// os << indent << "Process Mode: " << this->ProcessMode << "\n";
	os << indent << "Locator Max Level: " << this->LocatorMaxLevel << "\n";

	os << indent << "Capping: " << (this->Capping ? "On\n" : "Off\n");
	os << indent << "Cap Value: " << this->CapValue << "\n";
	// os << indent << "Process Mode: " << this->GetProcessModeAsString() << endl;
	os << indent << "Number Of Threads (for PerVoxel mode): " << this->NumberOfThreads << endl;
}


// //----------------------------------------------------------------------------
// // Templated append for VTK_VOXEL_MODE process mode and any type of output data
// template <class OT>
// void vtkSignedImplicitModellerAppendExecute(vtkSignedImplicitModeller *self,
// 											vtkDataSet *input, vtkImageData *outData,
// 											int outExt[6], double maxDistance, 
// 											vtkCellLocator *locator, int id, OT *)
// {
// 	int i, j, k;
// 	int subId;
// 	vtkIdType cellId;
// 	double pcoords[3];
// 	double *spacing, *origin;
// 	double maxDistance2 = maxDistance * maxDistance;   
// 	double x[3], prevDistance, prevDistance2, distance2, betterDistance;
// 	double closestPoint[3], mDist;

// 	// allocate weights for the EvaluatePosition
// 	double *weights = new double[input->GetMaxCellSize()];

// 	// Traverse each voxel; using CellLocator to find the closest point
// 	vtkGenericCell *cell = vtkGenericCell::New();

// 	spacing = outData->GetSpacing();
// 	origin = outData->GetOrigin();

// 	vtkImageProgressIterator<OT> outIt(outData, outExt, self, id);

// 	// so we know how to scale if desired
// 	double scaleFactor = 0; // 0 used to indicate not scaling
// 	double toDoubleScaleFactor = 0; // 0 used to indicate not scaling 
// 	double capValue = 0; // 0 used to indicate not clamping (float or double)
// 	if (self->GetOutputScalarType() != VTK_FLOAT && 
// 		self->GetOutputScalarType() != VTK_DOUBLE)
//     {
// 		capValue = self->GetCapValue();
// 		if (self->GetScaleToMaximumDistance())
// 		{
// 			scaleFactor = capValue / maxDistance;
// 			toDoubleScaleFactor = maxDistance / capValue;
// 		}
//     }

// 	int testIndex = 0;
// 	for (k = outExt[4]; k <= outExt[5]; k++) 
//     {
// 		x[2] = spacing[2] * k + origin[2];
// 		for (j = outExt[2]; j <= outExt[3]; j++)
// 		{
// 			cellId = -1;
// 			x[1] = spacing[1] * j + origin[1];
// 			OT* outSI = outIt.BeginSpan();
// 			for (i = outExt[0]; i <= outExt[1]; i++, testIndex++) 
// 			{
// 				x[0] = spacing[0] * i + origin[0];
// 				self->GetOutputDistance(*outSI, prevDistance, prevDistance2, toDoubleScaleFactor);
// 				betterDistance = -1;
        
// 				if (cellId != -1)
// 				{
// 					cell->EvaluatePosition(x, closestPoint, subId, pcoords,
// 										   distance2, weights);
// 					if (distance2 <= maxDistance2 && distance2 < prevDistance2)
// 					{
// 						mDist = sqrt(distance2);
// 						betterDistance = mDist;
// 					}
// 					else if (prevDistance2 < maxDistance2)
// 					{
// 						mDist = prevDistance;
// 					}
// 					else
// 					{
// 						mDist = maxDistance;
// 					}
// 				}
// 				else if (prevDistance2 < maxDistance2)
// 				{
// 					mDist = prevDistance;
// 				}
// 				else
// 				{
// 					mDist = maxDistance;
// 				}
        
// 				if (locator->FindClosestPointWithinRadius(x, mDist,
// 														  closestPoint, cell, cellId, subId, distance2) )
// 				{
// 					if(distance2 <= prevDistance2)
// 					{
// 						betterDistance = sqrt(distance2);
// 					}
// 				}
// 				else
// 				{
// 					cellId = -1;
// 				}

// 				if (betterDistance != -1)
// 				{
// 					self->SetOutputDistance(betterDistance, outSI, capValue, scaleFactor);
// 				}

// 				outSI++;
// 			}
// 			outIt.NextSpan();
// 		}
//     }
// 	cell->Delete();
// 	delete [] weights;
// }

// //----------------------------------------------------------------------------
// // This is the multithreaded piece of the append when doing per voxel
// // processing - it is called once for each thread, with each thread
// // taking a different slab of the output to work on.  The acutal work is done
// // in vtkSignedImplicitModellerAppendExecute; here we just setup for the per voxel
// // processing.
// static VTK_THREAD_RETURN_TYPE vtkSignedImplicitModeller_ThreadedAppend( void *arg )
// {
// 	int                      threadCount;
// 	int                      threadId;
// 	vtkSignedImplicitModellerAppendInfo *userData;
// 	vtkImageData *output;
// 	double maxDistance;
// 	int i;
// 	double *bounds, adjBounds[6];
// 	double *spacing;
// 	double *origin;
// 	int slabSize, slabMin, slabMax;
// 	int outExt[6];

// 	threadId = ((vtkMultiThreader::ThreadInfo *)(arg))->ThreadID;
// 	threadCount = ((vtkMultiThreader::ThreadInfo *)(arg))->NumberOfThreads;
// 	userData = (vtkSignedImplicitModellerAppendInfo *)
// 		(((vtkMultiThreader::ThreadInfo *)(arg))->UserData);

// 	if (userData->Input[threadId] == NULL)
//     {
// 		return VTK_THREAD_RETURN_VALUE;
//     }

// 	maxDistance = userData->MaximumDistance;

// 	output = userData->Modeller->GetOutput();
// 	spacing = output->GetSpacing();
// 	origin = output->GetOrigin();

// 	int *sampleDimensions = userData->Modeller->GetSampleDimensions();
// 	if (!output->GetPointData()->GetScalars())
//     {
// 		vtkGenericWarningMacro("Sanity check failed.");
// 		return VTK_THREAD_RETURN_VALUE;
//     }

// 	// break up into slabs based on threadId and threadCount
// 	slabSize = sampleDimensions[2] / threadCount;
// 	if (slabSize == 0) // in case threadCount >  sampleDimensions[2]
//     {
// 		slabSize = 1;
//     }
// 	slabMin = threadId * slabSize;
// 	if (slabMin >= sampleDimensions[2])
//     {
// 		return VTK_THREAD_RETURN_VALUE;
//     }
// 	slabMax = slabMin + slabSize - 1;
// 	if (threadId == threadCount - 1)
//     {
// 		slabMax = sampleDimensions[2] - 1;
//     }


// 	bounds = userData->Input[threadId]->GetBounds();
// 	for (i=0; i<3; i++)
//     {
// 		adjBounds[2*i] = bounds[2*i] - maxDistance;
// 		adjBounds[2*i+1] = bounds[2*i+1] + maxDistance;
//     }
  
// 	// compute dimensional bounds in data set
// 	for (i = 0; i < 3; i++)
//     {
// 		outExt[i*2] = (int) ((double)(adjBounds[2*i] - origin[i]) / 
// 							 spacing[i]);
// 		outExt[i*2+1] = (int) ((double)(adjBounds[2*i+1] - origin[i]) / 
// 							   spacing[i]);
// 		if (outExt[i*2] < 0)
// 		{
// 			outExt[i*2] = 0;
// 		}
// 		if (outExt[i*2+1] >= sampleDimensions[i])
// 		{
// 			outExt[i*2+1] = sampleDimensions[i] - 1;
// 		}
//     }

// 	// input not close enough to effect this slab
// 	if (outExt[4] > slabMax || outExt[5] < slabMin) 
//     {
// 		return VTK_THREAD_RETURN_VALUE;
//     }

// 	// adjust min/max to match slab
// 	if (outExt[4] < slabMin)
//     {
// 		outExt[4] = slabMin;
//     }
// 	if (outExt[5] > slabMax)
//     {
// 		outExt[5] = slabMax;
//     }
  
// 	vtkCellLocator *locator = vtkCellLocator::New();
  
// 	// Set up the cell locator.
// 	// If AutomaticOff, then NumberOfCellsPerBucket only used for allocating
// 	// memory.  If AutomaticOn, then NumberOfCellsPerBucket is used to guess
// 	// the depth for the uniform octree required to support
// 	// NumberOfCellsPerBucket (assuming uniform distribution of cells).
// 	locator->SetDataSet( userData->Input[threadId] );
// 	locator->AutomaticOff();
// 	locator->SetMaxLevel( userData->Modeller->GetLocatorMaxLevel() );
// 	locator->SetNumberOfCellsPerBucket( 1 );  
// 	locator->CacheCellBoundsOn();
// 	locator->BuildLocator();

// 	switch (userData->Modeller->GetOutputScalarType())
//     {
// 		vtkTemplateMacro(
// 			vtkSignedImplicitModellerAppendExecute( 
// 				userData->Modeller, 
// 				userData->Input[threadId], output, outExt, 
// 				userData->MaximumDistance, locator, threadId,
// 				static_cast<VTK_TT *>(0)));
//     default:
// 		vtkGenericWarningMacro("Execute: Unknown output ScalarType");
// 		return VTK_THREAD_RETURN_VALUE;
//     }

// 	locator->Delete();
// 	return VTK_THREAD_RETURN_VALUE;
// }

//----------------------------------------------------------------------------
// const char *vtkSignedImplicitModeller::GetProcessModeAsString()
// {
// 	if (this->ProcessMode == VTK_CELL_MODE)
//     {
// 		return "PerCell";
//     }
// 	else
//     {
// 		return "PerVoxel";
//     }
// }
