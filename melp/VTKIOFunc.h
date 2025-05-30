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

	template<int d>
	void OutputLagrangianParticlesAsVTPWithDisks(const LagrangeParticles<d>& particles, const fs::path& path) {
		Assert(d == 3, "OutputLagrangianParticlesAsVTUWithDisks is only implemented for 3D particles.");
		Info("Outputting Lagrangian particles as disk glyphs to VTU file: {}", path.string());

		const int N = particles.Size();

		// Points
		vtkNew<vtkPoints> points;
		points->SetNumberOfPoints(N);
		for (int i = 0; i < N; ++i) {
			const auto& x = particles.X(i);
			points->SetPoint(i, x[0], x[1], x[2]);
		}

		// Normals
		vtkNew<vtkFloatArray> normals;
		normals->SetNumberOfComponents(3);
		normals->SetNumberOfTuples(N);
		normals->SetName("Normals");
		for (int i = 0; i < N; ++i) {
			const auto& n = particles.Normal(i);
			normals->SetTuple3(i, n[0], n[1], n[2]);
		}




		// Scaling (radius based on area A, disk area = πr² => r = sqrt(A/π))
		vtkNew<vtkFloatArray> scale;
		scale->SetNumberOfComponents(1);
		scale->SetNumberOfTuples(N);
		scale->SetName("SA");
		for (int i = 0; i < N; ++i) {
			float r = std::sqrt(particles.SA(i) / pi);
			scale->SetValue(i, r);
		}

		// Set up source geometry: a unit disk
		vtkNew<vtkDiskSource> disk;
		disk->SetInnerRadius(0.0);     // solid disk
		disk->SetOuterRadius(1.0);     // will be scaled by Scale array
		disk->SetRadialResolution(1);  // low res for performance
		disk->SetCircumferentialResolution(20); // higher value for better visuals
		vtkNew<vtkTransform> transform;
		transform->RotateY(-90);
		vtkNew<vtkTransformPolyDataFilter> transformFilter;
		transformFilter->SetInputConnection(disk->GetOutputPort());
		transformFilter->SetTransform(transform);
		transformFilter->Update();

		// Set up glyphs
		vtkNew<vtkPolyData> glyph_input;
		glyph_input->SetPoints(points);
		// Set vectors explicitly
		glyph_input->GetPointData()->AddArray(normals);
		glyph_input->GetPointData()->SetVectors(normals);  // <- Add this!
		glyph_input->GetPointData()->SetNormals(normals);
		glyph_input->GetPointData()->AddArray(scale);

		//Set height
		vtkNew<vtkFloatArray> height;
		height->SetNumberOfComponents(1);
		height->SetNumberOfTuples(N);
		height->SetName("Height");
		for (int i = 0; i < N; ++i) {
			height->SetValue(i, particles.H(i));
		}
		glyph_input->GetPointData()->AddArray(height);

		vtkNew<vtkGlyph3D> glyphs;
		//glyphs->SetSourceConnection(disk->GetOutputPort());
		glyphs->SetSourceConnection(transformFilter->GetOutputPort());
		glyphs->SetInputData(glyph_input);
		glyphs->SetVectorModeToUseVector();       // <- Use "Normals" array
		glyphs->SetScaleModeToScaleByScalar();
		glyphs->SetColorModeToColorByScalar();
		glyphs->SetScaleFactor(1.0);
		glyphs->OrientOn();
		glyphs->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "SA");
		glyphs->Update();

		// Write to .vtu
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetFileName(path.string().c_str());
		writer->SetInputData(glyphs->GetOutput());
		writer->SetDataModeToBinary();
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