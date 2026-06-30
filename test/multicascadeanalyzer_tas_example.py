import numpy as np
import matplotlib.pyplot as plt

try:
    import pyvista as pv
except ImportError as exc:
    raise ImportError(
        "pyvista is required for 3D visualization. Install it with `pip install pyvista`."
    ) from exc

# Import components
from pyneutrace.components.nsource import NSource
from pyneutrace.components.spacer import Spacer
from pyneutrace.components.neutronguide import NeutronGuide
from pyneutrace.components.monochromator import Monochromator
from pyneutrace.components.vsample import VRodSample
from pyneutrace.components.multicascadeanalyzer import MultiCascadeAnalyzer
from pyneutrace.components.monitor import Monitor
from pyneutrace.components.virtualfilter import VirtualFilterAndMultiplier
from pyneutrace.instrument.instrumentassemble import InstrumentAssemble


def _add_pipeline_item(plotter, beamline, item, transform, fallback_color, show_edges):
    """Render one standard mesh-pipeline item in global coordinates."""
    item_type = item.get("type", "")
    transformed = beamline._transform_mesh_item(item, transform)

    if item_type == "structured_grid":
        grid = pv.StructuredGrid(
            transformed["data"]["X"],
            transformed["data"]["Y"],
            transformed["data"]["Z"],
        )
        plotter.add_mesh(
            grid,
            color=transformed.get("color", fallback_color),
            opacity=transformed.get("opacity", 0.5),
            show_edges=show_edges or transformed.get("show_edges", False),
        )
    elif item_type in ("polyline", "polyline_loop"):
        pts = np.column_stack((
            transformed["data"]["X"],
            transformed["data"]["Y"],
            transformed["data"]["Z"],
        ))
        if pts.shape[0] >= 2:
            line = pv.lines_from_points(pts, close=(item_type == "polyline_loop"))
            plotter.add_mesh(
                line,
                color=transformed.get("color", fallback_color),
                opacity=transformed.get("opacity", 0.8),
                line_width=transformed.get("line_width", 2),
            )
    elif item_type == "arrow":
        arrow = pv.Arrow(
            start=transformed["data"]["origin"],
            direction=transformed["data"]["vector"],
            scale="auto",
        )
        plotter.add_mesh(
            arrow,
            color=transformed.get("color", fallback_color),
            opacity=transformed.get("opacity", 1.0),
        )


def _translation_z(distance):
    transform = np.eye(4, dtype=float)
    transform[2, 3] = float(distance)
    return transform


def _build_analyzer_exit_monitor_entries(multicascade_entry):
    """Create visualization-only flat monitors at each analyzer exit plane."""
    multicascade = multicascade_entry["object"]
    multicascade_transform = multicascade_entry["transform"]
    monitor_entries = []

    for channel_index, (channel, channel_transform) in enumerate(
        zip(multicascade.channels, multicascade.channel_transforms)
    ):
        for analyzer_index, analyzer in enumerate(channel.analyzers):
            exit_transform = (
                multicascade_transform
                @ channel_transform
                @ _translation_z(analyzer.sample_R)
                @ analyzer.T_exit_from_entry
            )
            monitor_entries.append({
                "name": f"Analyzer Exit Monitor C{channel_index + 1} A{analyzer_index + 1}",
                "object": Monitor(
                    name=f"C{channel_index + 1}A{analyzer_index + 1}_ExitMonitor",
                    half_w=analyzer.exit_w / 2.0,
                    half_h=analyzer.exit_h / 2.0,
                ),
                "transform": exit_transform,
            })

    return monitor_entries


def plot_multicascade_exit_capture(multicascade):
    """Plot analyzer exit_capture data directly from stage_results."""
    n_stages = multicascade.num_analyzers_per_channel
    fig, axes = plt.subplots(1, n_stages, figsize=(6 * n_stages, 5), squeeze=False)
    axes = axes.ravel()

    for stage_idx, ax in enumerate(axes):
        x_all = []
        y_all = []
        w_all = []

        for channel_result in multicascade.channel_results:
            stage = channel_result["stage_results"][stage_idx]
            capture = stage["exit_capture"]
            if capture["x"].size == 0:
                continue
            x_all.append(capture["x"])
            y_all.append(capture["y"])
            w_all.append(capture["weight"])

        half_w = multicascade.exit_widths[stage_idx] / 2.0
        half_h = multicascade.exit_heights[stage_idx] / 2.0

        if x_all:
            x = np.concatenate(x_all)
            y = np.concatenate(y_all)
            w = np.concatenate(w_all)
            hist = ax.hist2d(
                x,
                y,
                bins=80,
                weights=w,
                range=[[-half_w, half_w], [-half_h, half_h]],
                cmap="inferno",
            )
            fig.colorbar(hist[3], ax=ax, label="Intensity")
            ax.set_title(
                f"Analyzer {stage_idx + 1} Exit Capture\n"
                f"count={x.size}, weight={w.sum():.3e}"
            )
        else:
            ax.text(0.5, 0.5, "No exit_capture", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"Analyzer {stage_idx + 1} Exit Capture")

        ax.set_xlim(-half_w, half_w)
        ax.set_ylim(-half_h, half_h)
        ax.set_xlabel("Exit x [mm]")
        ax.set_ylabel("Exit y [mm]")
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        f"{multicascade.name} exit_capture ({multicascade.exit_capture_mode} mode)",
        fontsize=12,
    )
    fig.tight_layout()
    plt.show()

def visualize_v4_beamline_3d(beamline, show_edges=False, extra_entries=None):
    """Render the full v4 multicascade beamline in global coordinates."""
    plotter = pv.Plotter(window_size=(1400, 900))
    plotter.add_text(f"3D Geometry Check: {beamline.name}", font_size=12)
    colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"]

    for idx, entry in enumerate(beamline.components):
        component = entry["object"]
        transform = entry["transform"]
        fallback_color = colors[idx % len(colors)]

        if hasattr(component, "get_mesh_pipeline"):
            for item in component.get_mesh_pipeline():
                _add_pipeline_item(plotter, beamline, item, transform, fallback_color, show_edges)
        elif type(component).__name__ == "VirtualFilterAndMultiplier":
            # This is a statistical helper stage, not a physical 3D object.
            # Skip it so the geometry plot shows only real instrument parts.
            continue
        else:
            raise NotImplementedError(
                f"Component '{entry['name']}' ({type(component).__name__}) does not provide 3D mesh data."
            )

    if extra_entries is not None:
        for entry in extra_entries:
            component = entry["object"]
            transform = entry["transform"]
            for item in component.get_mesh_pipeline():
                _add_pipeline_item(plotter, beamline, item, transform, "gold", show_edges)

    sample_entry = beamline.get_component("V-Rod Sample")
    sample = sample_entry["object"]
    sample_center_transform = sample_entry["transform"] @ sample.T_exit_from_entry
    sample_center = sample_center_transform[:3, 3]
    plotter.add_mesh(
        pv.Sphere(radius=20.0, center=sample_center),
        color="red",
        opacity=0.85,
        name="sample",
    )

    plotter.show_axes()
    plotter.show_grid()
    # Some IDE/terminal combinations return immediately from plotter.show(),
    # which makes the VTK window flash and then disappear when the script ends.
    # Keep the window alive explicitly in this test so the 3D layout can be inspected.
    plotter.show(auto_close=False)
    try:
        input("Press Enter to close the 3D window...")
    finally:
        plotter.close()


def run_assembled_simulation(n_neutrons=20000):
    print("=====================================================")
    print("Assembling and Visualizing v4 MultiCascade Beamline")
    print("=====================================================\n")

    mono_rotation = 15.0
    analyzer_rotations = [15.0, 20.0]
    center_scattering_angle = 70.0
    num_channels = 11
    channel_coverage_angle = 5.0
    dead_angle = 0.4

    beamline = InstrumentAssemble(name="Test Beamline Configuration v4 MultiCascade 3D")

    source = NSource(entry_w=100.0, entry_h=200.0, exit_w=100.0, exit_h=200.0, length=200.0)
    spacer1 = Spacer(entry_w=200.0, entry_h=200.0, exit_w=200.0, exit_h=200.0, length=400.0)
    guide = NeutronGuide(entry_w=200.0, entry_h=200.0, exit_w=200.0, exit_h=200.0, length=1000.0)
    spacer2 = Spacer(entry_w=200.0, entry_h=200.0, exit_w=200.0, exit_h=200.0, length=400.0)
    mono = Monochromator(
        entry_w=200.0,
        entry_h=200.0,
        exit_w=100.0,
        exit_h=200.0,
        length=800.0,
        rotation=mono_rotation,
    )
    spacer3 = Spacer(entry_w=100.0, entry_h=200.0, exit_w=100.0, exit_h=200.0, length=400.0)
    sample = VRodSample(entry_w=100.0, entry_h=200.0)
    virtual_filter = VirtualFilterAndMultiplier(
        min_lam=0.0,
        max_lam=10.0,
        enable_multiplier=True,
        multiplier=20,
        half_w=sample.entry_w / 2.0,
        half_h=sample.entry_h / 2.0,
        rng_seed=42,
    )
    multicascade = MultiCascadeAnalyzer(
        center_angle=center_scattering_angle,
        num_channels=num_channels,
        channel_coverage_angle=channel_coverage_angle,
        dead_angle=dead_angle,
        sample_R_list=[400.0, 700.0],
        analyzer_rotations=analyzer_rotations,
        entry_h=160.0,
        exit_w=220.0,
        exit_h=180.0,
        length=400.0,
        name="Multi-Cascade Analyzer",
    )
    post_sample_filter = VirtualFilterAndMultiplier(
        min_lam=0.0,
        max_lam=10.0,
        enable_multiplier=True,
        multiplier=10,
        half_w=multicascade.entry_w / 2.0,
        half_h=multicascade.entry_h / 2.0,
        rng_seed=84,
    )

    beamline.add_component("Neutron Source", source)
    beamline.add_component("Spacer 1", spacer1)
    beamline.add_component("Neutron Guide", guide)
    beamline.add_component("Spacer 2", spacer2)
    beamline.add_component("Monochromator", mono)
    beamline.add_component("Virtual Filter", virtual_filter)
    beamline.add_component("Spacer 3", spacer3)
    beamline.add_component("V-Rod Sample", sample)

    sample_entry = beamline.get_component("V-Rod Sample")
    sample_center_transform = sample_entry["transform"] @ sample.T_exit_from_entry
    beamline.add_component(
        "Post-Sample Virtual Filter",
        post_sample_filter,
        transform=sample_center_transform,
        auto_chain=False,
    )
    beamline.add_component("MultiCascade Analyzer", multicascade)

    print(f"Running demo instrument with {n_neutrons} rays...")
    print(f"Monochromator rotation: {mono_rotation:.1f} deg")
    print(f"Analyzer rotations: {analyzer_rotations}")
    print(f"MultiCascade center angle: {center_scattering_angle:.1f} deg")
    print(f"MultiCascade channels: {num_channels}")
    print(f"Channel coverage angle: {channel_coverage_angle:.1f} deg")
    print(f"Dead angle: {dead_angle:.1f} deg")
    print(f"First entry window [mm]: {multicascade.entry_w:.2f} x {multicascade.entry_h:.2f}")
    print(f"Entry widths [mm]: {np.round(multicascade.entry_widths, 3)}")
    print(f"Exit widths [mm]: {np.round(multicascade.exit_widths, 3)}")
    print(f"Exit heights [mm]: {np.round(multicascade.exit_heights, 3)}")
    print(f"Analyzer radii [mm]: {np.round(multicascade.sample_R_list, 3)}")

    pos, vel, t, w = source.generate_beam(n_neutrons)
    speed = np.linalg.norm(vel, axis=1)

    transport_to_mono = spacer1.length + guide.length + spacer2.length
    rng = np.random.default_rng(42)
    target_x = rng.uniform(-mono.entry_w / 2.0, mono.entry_w / 2.0, n_neutrons)
    target_y = rng.uniform(-mono.entry_h / 2.0, mono.entry_h / 2.0, n_neutrons)
    direction_to_mono = np.column_stack((
        target_x - pos[:, 0],
        target_y - pos[:, 1],
        np.full(n_neutrons, transport_to_mono),
    ))
    direction_to_mono /= np.linalg.norm(direction_to_mono, axis=1, keepdims=True)
    vel = direction_to_mono * speed[:, np.newaxis]

    final_pos, final_vel, final_t, final_w = beamline.run(
        initPos=pos,
        initVel=vel,
        initTime=t,
        initWeight=w,
    )

    beamline.print_summary()
    print(f"\nResult: {len(final_pos)} rays reached the end of the beamline.")
    print(f"Total weight after propagation: {np.sum(final_w):.4e}")
    for channel_segment in multicascade.output_segments:
        ch_idx = channel_segment["channel_index"]
        channel_result = multicascade.channel_results[ch_idx]
        print(
            f"Channel {ch_idx} @ {channel_segment['channel_angle']:.2f} deg: "
            f"{channel_segment['count']} reflected neutrons, "
            f"total weight {channel_result['reflected_weight']:.4e}"
        )
        for stage in channel_result["stage_results"]:
            print(
                f"  Analyzer {stage['analyzer_index'] + 1}: "
                f"input={stage['input_count']}, "
                f"reflected={stage['reflected_count']}, "
                f"weight={stage['reflected_weight']:.4e}"
            )

    analyzer_exit_monitors = _build_analyzer_exit_monitor_entries(
        beamline.get_component("MultiCascade Analyzer")
    )
    print(f"Analyzer exit monitors added for 3D view: {len(analyzer_exit_monitors)}")
    plot_multicascade_exit_capture(multicascade)

    visualize_v4_beamline_3d(
        beamline,
        show_edges=False,
        extra_entries=analyzer_exit_monitors,
    )


if __name__ == "__main__":
    run_assembled_simulation()
