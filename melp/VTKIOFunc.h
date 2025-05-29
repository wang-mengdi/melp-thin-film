#pragma once

#include "LagrangeParticles.h"
#include "EulerParticles.h"

#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkDiskSource.h>

namespace fs = std::filesystem;

namespace VTKIOFunc {

	template<int d>
	void OutputLagrangianParticlesAsVTU(const LagrangeParticles<d>& particles, const fs::path& path) {
		Assert(d == 3, "OutputLagrangianParticlesAsVTU is only implemented for 3D particles.");
		Info("Outputting Lagrangian particles to VTU file: {}\n", path.string());

		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		// Add points
		vtkNew<vtkFloatArray> positions;
		positions->SetNumberOfComponents(3);  // 3D points
		positions->SetNumberOfTuples(particles.Size());
		positions->SetName("Positions");
		for (int i = 0; i < particles.Size(); ++i) {
			const auto pos = particles.X(i);
			positions->SetTuple3(i, pos[0], pos[1], pos[2]);
		}
		vtkNew<vtkPoints> nodes;
		nodes->SetData(positions);
		unstructured_grid->SetPoints(nodes);

		// Add height
		vtkNew<vtkFloatArray> heights;
		heights->SetNumberOfComponents(1);
		heights->SetNumberOfTuples(particles.Size());
		heights->SetName("Height");
		for (int i = 0; i < particles.Size(); ++i) {
			heights->SetValue(i, particles.H(i));
		}
		unstructured_grid->GetPointData()->AddArray(heights);

		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		// Write the output file
		writer->SetFileName(path.string().c_str());
		writer->SetInputData(unstructured_grid);
		writer->SetDataModeToBinary();  // Optional: Use binary mode for smaller file size
		writer->Write();
	}

	template<int d>
	void OutputEulerianParticlesAsVTU(const EulerParticles<d>& particles, const fs::path& path) {
		Assert(d == 3, "OutputEulerianParticlesAsVTU is only implemented for 3D particles.");
		Info("Outputting Eulerian particles to VTU file: {}\n", path.string());

		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		// Add points
		vtkNew<vtkFloatArray> positions;
		positions->SetNumberOfComponents(3);  // 3D points
		positions->SetNumberOfTuples(particles.Size());
		positions->SetName("Positions");
		for (int i = 0; i < particles.Size(); ++i) {
			const auto pos = particles.X(i);
			positions->SetTuple3(i, pos[0], pos[1], pos[2]);
		}
		vtkNew<vtkPoints> nodes;
		nodes->SetData(positions);
		unstructured_grid->SetPoints(nodes);

		// Add height
		vtkNew<vtkFloatArray> heights;
		heights->SetNumberOfComponents(1);
		heights->SetNumberOfTuples(particles.Size());
		heights->SetName("Height");
		for (int i = 0; i < particles.Size(); ++i) {
			heights->SetValue(i, particles.H(i));
		}
		unstructured_grid->GetPointData()->AddArray(heights);

		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		// Write the output file
		writer->SetFileName(path.string().c_str());
		writer->SetInputData(unstructured_grid);
		writer->SetDataModeToBinary();  // Optional: Use binary mode for smaller file size
		writer->Write();
	}

}