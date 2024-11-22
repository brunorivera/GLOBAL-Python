from clases import Detector, Radiacion, Virus, Sanador
from typing import List


def mostrar_matriz(matriz: List[List[str]]) -> None:
    """
    Muestra la matriz de ADN en formato tabular, separada de otros elementos.
    """
    if not matriz:
        print("La matriz está vacía.")
        return

    print("\n" + "=" * 30)
    print("MATRIZ DE ADN:")
    print("=" * 30)
    print("    " + "  ".join(str(i) for i in range(len(matriz[0]))))
    for idx, fila in enumerate(matriz):
        print(f"{idx} | " + "  ".join(fila))
    print("=" * 30 + "\n")


def ingresar_matriz() -> List[List[str]]:
    """
    Solicita al usuario que ingrese una matriz de ADN.
    """
    while True:
        try:
            print("\nIngrese una matriz de ADN (6 filas de 6 bases nitrogenadas, separadas por una coma y un espacio):")
            print("Ejemplo: AGATCA, GATTCA, CAACAT, GAGCTA, ATTGCG, CTGTTC")
            entrada = input("Matriz: ").strip().upper()
            elementos = entrada.split(", ")
            if len(elementos) != 6 or not all(len(fila) == 6 for fila in elementos):
                raise ValueError("Formato incorrecto. Debe ingresar 6 filas de 6 bases nitrogenadas.")
            if not all(base in "ATCG" for fila in elementos for base in fila):
                raise ValueError("Bases no válidas. Solo se permiten A, T, C y G.")
            return [list(fila) for fila in elementos]
        except ValueError as e:
            print(f"Error: {e}")


def pedir_entero(mensaje: str, min_val: int, max_val: int) -> int:
    """
    Solicita un número entero dentro de un rango válido.
    """
    while True:
        try:
            valor = int(input(mensaje))
            if min_val <= valor <= max_val:
                return valor
            else:
                print(f"Error: El número debe estar entre {min_val} y {max_val}.")
        except ValueError:
            print("Error: Debe ingresar un número entero válido.")


def validar_base_nitrogenada() -> str:
    """
    Valida que la base nitrogenada ingresada sea A, T, C o G.
    """
    while True:
        base = input("Base nitrogenada (A/T/C/G): ").strip().upper()
        if base and base in "ATCG":
            return base
        elif not base:
            print("Error: El campo está vacío. Debe ingresar una base nitrogenada (A, T, C o G).")
        else:
            print("Error: Base inválida. Solo se permiten A, T, C y G. Intente nuevamente.")


def validar_tipo_mutacion() -> str:
    """
    Valida que el tipo de mutación sea H, V o D.
    """
    while True:
        tipo = input("Tipo de mutación (H para Horizontal, V para Vertical, D para Diagonal): ").strip().upper()
        if tipo and tipo in "HVD":
            return tipo
        elif not tipo:
            print("Error: El campo está vacío. Debe ingresar un tipo de mutación (H, V o D).")
        else:
            print("Error: Tipo inválido. Solo se permiten H, V o D. Intente nuevamente.")


def imprimir_resultado(booleano: bool, mensaje: str) -> None:
    """
    Imprime el resultado booleano y el mensaje asociado.
    """
    print(f"\nResultado booleano: {booleano}")
    print(mensaje)
    
#Main

def main() -> None:
    """
    Función principal del programa.
    """
    matriz = ingresar_matriz()
    mostrar_matriz(matriz)

    opciones = {
        "1": "Detectar mutaciones",
        "2": "Crear mutantes",
        "3": "Sanar ADN",
        "4": "Ingresar una nueva matriz",
        "5": "Salir",
    }

    while True:
        print("\nOpciones:")
        for key, val in opciones.items():
            print(f"{key}. {val}")
        opcion = input("\nSeleccione una opción: ").strip()

        if opcion == "1":
            detector = Detector(matriz)
            resultado = detector.detectar_mutantes()
            print(f"\nResultado booleano: {resultado}")
            if resultado:
                print("Mutaciones detectadas en las siguientes posiciones:")
                for fila, columna, tipo in detector.posiciones_mutaciones:
                    print(f"- Fila: {fila}, Columna: {columna}, Tipo: {tipo}")
            else:
                print("No se detectaron mutaciones.")

        elif opcion == "2":
            base = validar_base_nitrogenada()
            while True:
                fila = pedir_entero("Fila inicial (0-5): ", 0, 5)
                columna = pedir_entero("Columna inicial (0-5): ", 0, 5)
                tipo = validar_tipo_mutacion()
                try:
                    if tipo in ["H", "V"]:
                        radiacion = Radiacion(base, 4, tipo)
                        matriz = radiacion.crear_mutante(matriz, (fila, columna), tipo)
                    elif tipo == "D":
                        virus = Virus(base, 4, "Diagonal")
                        matriz = virus.crear_mutante(matriz, (fila, columna))
                    print("\nMATRIZ MUTADA:")
                    break
                except ValueError as e:
                    print(f"Error: {e}. Intente nuevamente.")

        elif opcion == "3":
            sanador = Sanador(matriz)
            matriz_sanada, cambios_realizados = sanador.sanar_mutantes(matriz)
            if cambios_realizados:
                print("\nLa matriz ha sido sanada correctamente.")
            else:
                print("\nLa matriz ya estaba sana. No se realizaron cambios.")
            matriz = matriz_sanada

        elif opcion == "4":
            matriz = ingresar_matriz()
            print("\nNueva matriz ingresada:")
            mostrar_matriz(matriz)

        elif opcion == "5":
            print("\nSaliendo...")
            break

        else:
            print("Opción inválida. Intente nuevamente.")

        if opcion != "4":  # Evitar duplicar la matriz tras ingreso
            mostrar_matriz(matriz)


if __name__ == "__main__":
    main()
