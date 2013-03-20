/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSignedImplicitModeller.h

  Copyright (c) Alex Mansfield
  All rights reserved.
  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
// .NAME vtkSignedImplicitModeller - compute signed distance from input geometry on structured point dataset
// .SECTION Description
// vtkSignedImplicitModeller is a filter ...

#ifndef __vtkSignedImplicitModeller_h
#define __vtkSignedImplicitModeller_h

#include <list>
#include <vector>
#include <math.h>

#include "vtkFiltersHybridModule.h" // For export macro
#include "vtkImageAlgorithm.h"
#include "vtkMath.h"
#include "vtkPolyData.h"

/* #define VTK_VOXEL_MODE   0 */
/* #define VTK_CELL_MODE    1 */

class vtkDataArray;
class vtkExtractGeometry;
class vtkMultiThreader;

class VTKFILTERSHYBRID_EXPORT vtkSignedImplicitModeller : public vtkImageAlgorithm 
{
public:
	class Barycentric
	{
	protected:
		bool initialized;
		double origin[3], dir0[3], dir1[3], p0_plane[2], p1_plane[2], p2_plane[2];
		double det;

		static void normalize(double * const x)
		{
			double length = 0.0;
			for(int i=0; i<3; i++)
				length += x[i]*x[i];

			length = std::sqrt(length);
			for(int i=0; i<3; i++)
				x[i] /= length;
			
			return;
		}

		void planar(double const * const x, double * const p) const
		{
			for(int i=0; i<2; i++)
				p[i] = 0.0;

			for(int i=0; i<3; i++)
			{
				p[0] += dir0[i] *(x[i] - origin[i]);
				p[1] += dir1[i] *(x[i] - origin[i]);
			}

			return;
		}

	public:

		Barycentric(double const * const p0, double const * const p1, double const * const p2)
		{
			// compute origin
			for(int i=0; i<3; i++)
				origin[i] = p2[i];
			
			// compute first direction in the plane, dir0
			for(int i=0; i<3; i++)
				dir0[i] = p1[i] - origin[i];
			normalize(dir0);

			// compute second direction in the plane, dir1
			double dot_product = 0.0;
			for(int i=0; i<3; i++)
			{
				dir1[i] = p0[i] - origin[i];
				dot_product += dir0[i]*dir1[i];
			}
			for(int i=0; i<3; i++)
				dir1[i] -= dot_product*dir0[i];
			normalize(dir1);

			// p2_plane is always (0,0)
			p2_plane[0] = 0.0;
			p2_plane[1] = 0.0;
			// p1_plane is always (a,0)
			planar(p1, p1_plane);
			planar(p0, p0_plane);

			// compute det
			det = -1.0*p0_plane[1]*p1_plane[0];

			// check that matrix is invertible
			if( (det*det) < 1e-9 )
				initialized = false;
			else
				initialized = true;
		}

		void compute(double const * const x, double * const b) const
	   	{
			for(int i=0; i<3; i++)
				b[i] = 0.0;

			if(initialized)
			{
				// convert x to planar co-ordinates
				double x_plane[2];
				planar(x, x_plane);

				// solve for Barycentric co-ordinates
				b[0] = vtkMath::Determinant2x2(x_plane  , p1_plane) / det;
				b[1] = vtkMath::Determinant2x2(p0_plane , x_plane)  / det;
				b[2] = 1.0 - (b[0] + b[1]);
			}

			return;
		}

		bool get_initialized() const
		{
			return initialized;
		}
	};

	vtkTypeMacro(vtkSignedImplicitModeller,vtkImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// Construct with sample dimensions=(50,50,50), and so that model bounds are
	// automatically computed from the input. Capping is turned on with CapValue
	// equal to a large positive number.
	static vtkSignedImplicitModeller *New();

	// Description:
	// Compute ModelBounds from input geometry. If input is not specified, the
	// input of the filter will be used.
	double ComputeModelBounds(vtkDataSet *input = NULL);

	// Description:
	// Set/Get the voxel spacing.
	vtkSetMacro(Spacing,double);
	vtkGetMacro(Spacing,double);

	// Description:
	// Set/Get the i-j-k dimensions on which to sample distance function.
	vtkGetVectorMacro(SampleDimensions,int,3);
	void SetSampleDimensions(int i, int j, int k);
	void SetSampleDimensions(int dim[3]);

	// Description:
	// Set / get the distance away from surface of input geometry to
	// sample. Smaller values make large increases in performance.
	vtkSetClampMacro(MaximumDistance,double,0.0,1.0);
	vtkGetMacro(MaximumDistance,double);

	// Description:
	// Set / get the region in space in which to perform the sampling. If
	// not specified, it will be computed automatically.
	vtkSetVector6Macro(ModelBounds,double);
	vtkGetVectorMacro(ModelBounds,double,6);

	// Description:
	// Control how the model bounds are computed. If the ivar AdjustBounds
	// is set, then the bounds specified (or computed automatically) are modified
	// by the fraction given by AdjustDistance. This means that the model
	// bounds is expanded in each of the x-y-z directions.
	vtkSetMacro(AdjustBounds,int);
	vtkGetMacro(AdjustBounds,int);
	vtkBooleanMacro(AdjustBounds,int);
  
	// Description:
	// Specify the amount to grow the model bounds (if the ivar AdjustBounds
	// is set). The value is a fraction of the maximum length of the sides
	// of the box specified by the model bounds.
	vtkSetClampMacro(AdjustDistance,double,-1.0,1.0);
	vtkGetMacro(AdjustDistance,double);

	// Description:
	// Control how the model bounds are computed. If the ivar AdjustBoundsAbsolute
	// is set, then the bounds specified (or computed automatically) are modified
	// by the absolute amount given by AdjustDistanceAbsolute. This means that the model
	// bounds is expanded in each of the x-y-z directions.
	vtkSetMacro(AdjustBoundsAbsolute,int);
	vtkGetMacro(AdjustBoundsAbsolute,int);
	vtkBooleanMacro(AdjustBoundsAbsolute,int);
  
	// Description:
	// Specify the amount to grow the model bounds (if the ivar AdjustBoundsAbsolute
	// is set). Applied after the relative adjustment described by AdjustDistance.
	vtkSetMacro(AdjustDistanceAbsolute,double);
	vtkGetMacro(AdjustDistanceAbsolute,double);

	// Description:
	// Specify whether to compute SDF from complete or segmented mesh. SDF only computed using
	// points on the segmented with color (255,255,255).
	void SetBoxesSegmentation(std::list<std::vector<double> > const & _BoxesSegmentation);
	void GetBoxesSegmentation(std::list<std::vector<double> > & _BoxesSegmentation);

	// Description:
	// The outer boundary of the structured point set can be assigned a 
	// particular value. This can be used to close or "cap" all surfaces.
	vtkSetMacro(Capping,int);
	vtkGetMacro(Capping,int);
	vtkBooleanMacro(Capping,int);
  
	// Description:
	// Specify the capping value to use. The CapValue is also used as an
	// initial distance value at each point in the dataset.
	void SetCapValue(double value);
	vtkGetMacro(CapValue,double);

	// Description:
	// If a non-floating output type is specified, the output distances can be
	// scaled to use the entire positive scalar range of the output type 
	// specified (up to the CapValue which is equal to the max for the type 
	// unless modified by the user).  For example, if ScaleToMaximumDistance
	// is On and the OutputScalarType is UnsignedChar the distances saved in the
	// output would be linearly scaled between 0 (for distances "very close" to
	// the surface) and 255 (at the specifed maximum distance)... assuming the 
	// CapValue is not changed from 255.
	vtkSetMacro(ScaleToMaximumDistance, int);
	vtkGetMacro(ScaleToMaximumDistance, int);
	vtkBooleanMacro(ScaleToMaximumDistance,int);

	// Description:
	// Specify whether to visit each cell once per append or each voxel once
	// per append.  Some tests have shown once per voxel to be faster
	// when there are a lot of cells (at least a thousand?); relative
	// performance improvement increases with addition cells.  Primitives
	// should not be stripped for best performance of the voxel mode.  
	/* vtkSetClampMacro(ProcessMode, int, 0, 1); */
	/* vtkGetMacro(ProcessMode, int); */
	/* void SetProcessModeToPerVoxel() {this->SetProcessMode(VTK_VOXEL_MODE);} */
	/* void SetProcessModeToPerCell()  {this->SetProcessMode(VTK_CELL_MODE);} */
	/* const char *GetProcessModeAsString(void); */

	// Description:
	// Specify the level of the locator to use when using the per voxel
	// process mode.
	vtkSetMacro(LocatorMaxLevel,int);
	vtkGetMacro(LocatorMaxLevel,int);

	// Description:
	// Set / Get the number of threads used during Per-Voxel processing mode
	vtkSetClampMacro( NumberOfThreads, int, 1, VTK_MAX_THREADS );
	vtkGetMacro( NumberOfThreads, int );

	// Description:
	// Set the desired output scalar type.
	void SetOutputScalarType(int type);
	vtkGetMacro(OutputScalarType,int);
	void SetOutputScalarTypeToFloat(){this->SetOutputScalarType(VTK_FLOAT);};
	void SetOutputScalarTypeToDouble(){this->SetOutputScalarType(VTK_DOUBLE);};
	void SetOutputScalarTypeToInt(){this->SetOutputScalarType(VTK_INT);};
	void SetOutputScalarTypeToUnsignedInt()
    {this->SetOutputScalarType(VTK_UNSIGNED_INT);};
	void SetOutputScalarTypeToLong(){this->SetOutputScalarType(VTK_LONG);};
	void SetOutputScalarTypeToUnsignedLong()
    {this->SetOutputScalarType(VTK_UNSIGNED_LONG);};
	void SetOutputScalarTypeToShort(){this->SetOutputScalarType(VTK_SHORT);};
	void SetOutputScalarTypeToUnsignedShort()   
    {this->SetOutputScalarType(VTK_UNSIGNED_SHORT);};
	void SetOutputScalarTypeToUnsignedChar()
    {this->SetOutputScalarType(VTK_UNSIGNED_CHAR);};
	void SetOutputScalarTypeToChar()
    {this->SetOutputScalarType(VTK_CHAR);};

	// Description:
	// Initialize the filter for appending data. You must invoke the
	// StartAppend() method before doing successive Appends(). It's also a
	// good idea to manually specify the model bounds; otherwise the input
	// bounds for the data will be used.
	void StartAppend();

	// Description:
	// Append a data set to the existing output. To use this function,
	// you'll have to invoke the StartAppend() method before doing
	// successive appends. It's also a good idea to specify the model
	// bounds; otherwise the input model bounds is used. When you've
	// finished appending, use the EndAppend() method.
	void Append(vtkDataSet *input);

	// Description:
	// Method completes the append process.
	void EndAppend();

	// See the vtkAlgorithm for a desciption of what these do
	int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

	template<class OT>
	void SetOutputDistance(double distance, OT *outputValue, double capValue=0, double scaleFactor=0);
	template<class OT>
	void GetOutputDistance(OT outputValue, double & distance, double & distance2, double scaleFactor=0);

	// computation of signed distance function
	int ComputeSignedDistance(vtkPolyData * input, int cellNum, Barycentric const & barycentric, double x[3], double & sdist);

protected:
	vtkSignedImplicitModeller();
	~vtkSignedImplicitModeller();

	double GetScalarTypeMax(int type);

	virtual int RequestInformation (vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestData (vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
	void StartAppend(int internal);
	void Cap(vtkDataArray *s);

	vtkMultiThreader *Threader;
	int              NumberOfThreads;

	double Spacing;
	int SampleDimensions[3];
	double MaximumDistance;
	double ModelBounds[6];
	int Capping;
	double CapValue;
	int DataAppended;
	int AdjustBounds, AdjustBoundsAbsolute;
	double AdjustDistance, AdjustDistanceAbsolute;
	std::list<std::vector<double> > BoxesSegmentation;

	/* int ProcessMode; */
	int LocatorMaxLevel;
	int OutputScalarType;
	int ScaleToMaximumDistance;

	// flag to limit to one ComputeModelBounds per StartAppend
	int BoundsComputed; 
  
	// the max distance computed during that one call
	double InternalMaxDistance; 

	virtual int FillInputPortInformation(int, vtkInformation*);
	
private:
	vtkSignedImplicitModeller(const vtkSignedImplicitModeller&);  // Not implemented.
	void operator=(const vtkSignedImplicitModeller&);  // Not implemented.
};

#endif


