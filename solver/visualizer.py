import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def read_csv(filename):
    rectangles = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, delimiter=";")
        for row in reader:
            rect = {
                'x': int(row['x']),
                'y': int(row['y']),
                'width': int(row['width']),
                'height': int(row['height'])
            }
            rectangles.append(rect)
    return rectangles

def plot_rectangles(rectangles):
    fig, ax = plt.subplots()
    for rect in rectangles:
        ax.add_patch(patches.Rectangle(
            (rect['x'], rect['y']), rect['width'], rect['height'],
            edgecolor='blue', facecolor='none', lw=2
        ))
    ax.set_aspect('equal', 'box')
    ax.set_xlim((-1, 21))
    ax.set_ylim((-1, 31))
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Rectangles from CSV')
    plt.show()

if __name__ == "__main__":
    filename = './output.csv'
    rectangles = read_csv(filename)
    plot_rectangles(rectangles)
