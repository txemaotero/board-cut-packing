import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def read_csv(filename):
    rectangles = []
    with open(filename, "r") as file:
        reader = csv.DictReader(file, delimiter=";")
        for row in reader:
            rect = {
                "x": int(row["x"]),
                "y": int(row["y"]),
                "width": int(row["width"]),
                "height": int(row["height"]),
            }
            rectangles.append(rect)
    return rectangles


def plot_rectangles(rectangles, title):
    fig, ax = plt.subplots()

    # Colores predefinidos
    colors = [
        "red",
        "green",
        "blue",
        "orange",
        "purple",
        "brown",
        "pink",
        "gray",
        "yellow",
        "cyan",
    ]

    if not rectangles:
        raise ValueError("La lista de rectángulos no puede estar vacía")

    # Dibuja el primer rectángulo con borde negro y sin relleno
    first_rect = rectangles[0]
    ax.add_patch(
        patches.Rectangle(
            (first_rect["x"], first_rect["y"]),
            first_rect["width"],
            first_rect["height"],
            edgecolor="black",
            facecolor="none",
            ls="--",
            lw=1,
        )
    )

    # Dibuja los demás rectángulos con colores diferentes y relleno transparente
    for i, rect in enumerate(rectangles[1:], start=1):
        color = colors[i % len(colors)]  # Usar un color de la lista de colores
        ax.add_patch(
            patches.Rectangle(
                (rect["x"], rect["y"]),
                rect["width"],
                rect["height"],
                edgecolor="none",
                facecolor=color,
                alpha=0.5,
                lw=2,
            )
        )
        rotation = "horizontal"
        if rect["height"] > rect["width"]:
            rotation = "vertical"
        ax.text(
            rect["x"] + rect["width"] / 2,
            rect["y"] + rect["height"] / 2,
            f"{rect['width']} x {rect['height']}",
            ha="center",
            va="center",
            rotation=rotation,
            fontsize=10,
            color="black",
        )

    # Calcula los límites basándose en el primer rectángulo
    margin = 100  # margen adicional para la visualización
    xlim = (first_rect["x"] - margin, first_rect["x"] + first_rect["width"] + margin)
    ylim = (first_rect["y"] - margin, first_rect["y"] + first_rect["height"] + margin)

    ax.set_aspect("equal", "box")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    ax.axis('off')
    title = title.replace("board", "tablero_").replace(".csv", "")
    plt.title(title)
    plt.savefig(title + ".pdf")
    plt.close(fig)


if __name__ == "__main__":
    from glob import glob

    filenames = glob("outputs/board*mm_*.csv")
    # filenames = ["./problematic.csv"]
    for fname in filenames:
        rectangles = read_csv(fname)
        plot_rectangles(rectangles, fname)

