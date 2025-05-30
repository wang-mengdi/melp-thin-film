#pragma once

#include "LagrangeParticles.h"
#include "EulerParticles.h"
#include "Constants.h"

#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkDiskSource.h>
#include <vtkGlyph3D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

namespace fs = std::filesystem;

namespace VTKIOFunc {

	void OutputCircleDisksAsVTP(const Array<Vector<real, 3>>& point_array, const Array<Vector<real, 3>>& normal_array, const Array<real>& area_array, const Array<real>& height_array, const fs::path& path);
	void OutputPointCloudDataAsVTP(const Array<Vector3>& positions, const std::vector<std::pair<Array<real>, std::string>> scalar_fields, const std::vector<std::pair<Array<Vector3>, std::string>> vector_fields, const fs::path& path);
}