from typing import List, Tuple
import random


class Detector:
    """
    Clase para detectar mutaciones en una matriz de ADN.
    """
    def __init__(self, matriz: List[List[str]]):
        """
        Inicializa un detector con la matriz de ADN proporcionada.
        """
        self.matriz = matriz
        self.mutaciones_encontradas = 0
        self.posiciones_mutaciones = []

    def detectar_mutantes(self, matriz: List[List[str]] = None) -> bool:
        """
        Detecta mutaciones horizontales, verticales y diagonales en la matriz de ADN.
        Si encuentra al menos una mutación, devuelve True, de lo contrario False.

        Args:
            matriz (List[List[str]]): Matriz de ADN a analizar. Si no se proporciona, se utiliza la matriz del objeto.

        Returns:
            bool: True si se detectan mutaciones, False en caso contrario.
        """
        matriz = matriz or self.matriz
        self.mutaciones_encontradas = 0
        self.posiciones_mutaciones = []
        self._detectar_horizontal(matriz)
        self._detectar_vertical(matriz)
        self._detectar_diagonal(matriz)
        return self.mutaciones_encontradas > 0

    def _detectar_horizontal(self, matriz: List[List[str]]) -> None:
        """
        Detecta mutaciones horizontales en la matriz de ADN.
        """
        for fila in range(len(matriz)):
            for columna in range(len(matriz[fila]) - 3):
                secuencia = matriz[fila][columna:columna + 4]
                if len(set(secuencia)) == 1:
                    self.mutaciones_encontradas += 1
                    self.posiciones_mutaciones.append((fila, columna, "Horizontal"))

    def _detectar_vertical(self, matriz: List[List[str]]) -> None:
        """
        Detecta mutaciones verticales en la matriz de ADN.
        """
        for columna in range(len(matriz[0])):
            for fila in range(len(matriz) - 3):
                secuencia = [matriz[fila + i][columna] for i in range(4)]
                if len(set(secuencia)) == 1:
                    self.mutaciones_encontradas += 1
                    self.posiciones_mutaciones.append((fila, columna, "Vertical"))

    def _detectar_diagonal(self, matriz: List[List[str]]) -> None:
        """
        Detecta mutaciones diagonales en la matriz de ADN.
        """
        n = len(matriz)
        for diag in range(-n + 1, n):
            diag_principal = []
            diag_secundaria = []
            for i in range(max(0, diag), min(n, n + diag)):
                diag_principal.append(matriz[i][i - diag])
                diag_secundaria.append(matriz[i][n - 1 - (i - diag)])
            for i in range(len(diag_principal) - 3):
                secuencia_p = diag_principal[i:i + 4]
                secuencia_s = diag_secundaria[i:i + 4]
                if len(set(secuencia_p)) == 1:
                    self.mutaciones_encontradas += 1
                    self.posiciones_mutaciones.append((i, i - diag, "Diagonal Principal"))
                if len(set(secuencia_s)) == 1:
                    self.mutaciones_encontradas += 1
                    self.posiciones_mutaciones.append((i, n - 1 - (i - diag), "Diagonal Secundaria"))


class Mutador:
    """
    Clase base abstracta para mutaciones de ADN.
    """
    def __init__(self, base_nitrogenada: str, longitud_mutacion: int, tipo_mutacion: str):
        """
        Inicializa el mutador con los atributos base nitrogenada, longitud y tipo.
        """
        self.base_nitrogenada = base_nitrogenada
        self.longitud_mutacion = longitud_mutacion
        self.tipo_mutacion = tipo_mutacion

    def crear_mutante(self):
        raise NotImplementedError("Debe implementarse en la clase derivada")


class Radiacion(Mutador):
    """
    Clase para mutaciones de ADN por radiación (horizontal o vertical).
    """
    def crear_mutante(self, matriz: List[List[str]], posicion: Tuple[int, int], orientacion: str) -> List[List[str]]:
        """
        Aplica una mutación de radiación en la matriz de ADN.

        Args:
            matriz (List[List[str]]): Matriz de ADN.
            posicion (Tuple[int, int]): Posición inicial de la mutación.
            orientacion (str): Orientación de la mutación ("H" o "V").

        Returns:
            List[List[str]]: Matriz de ADN mutada.
        """
        fila, columna = posicion
        if orientacion == "H" and columna + self.longitud_mutacion <= len(matriz[fila]):
            for i in range(self.longitud_mutacion):
                matriz[fila][columna + i] = self.base_nitrogenada
        elif orientacion == "V" and fila + self.longitud_mutacion <= len(matriz):
            for i in range(self.longitud_mutacion):
                matriz[fila + i][columna] = self.base_nitrogenada
        else:
            raise ValueError("La mutación se desbordaría de la matriz.")
        return matriz


class Virus(Mutador):
    """
    Clase para mutaciones de ADN por virus (diagonal).
    """
    def crear_mutante(self, matriz: List[List[str]], posicion: Tuple[int, int]) -> List[List[str]]:
        """
        Aplica una mutación de virus en la matriz de ADN.

        Args:
            matriz (List[List[str]]): Matriz de ADN.
            posicion (Tuple[int, int]): Posición inicial de la mutación.

        Returns:
            List[List[str]]: Matriz de ADN mutada.
        """
        fila, columna = posicion
        if fila + self.longitud_mutacion <= len(matriz) and columna + self.longitud_mutacion <= len(matriz[0]):
            for i in range(self.longitud_mutacion):
                matriz[fila + i][columna + i] = self.base_nitrogenada
        else:
            raise ValueError("La mutación se desbordaría de la matriz.")
        return matriz


class Sanador:
    """
    Clase para sanar mutaciones en una matriz de ADN.
    """
    def __init__(self, matriz: List[List[str]]):
        """
        Inicializa el sanador con la matriz original.
        """
        self.matriz_original = matriz

    def sanar_mutantes(self, matriz: List[List[str]]) -> Tuple[List[List[str]], bool]:
        """
        Sana una matriz de ADN generando una nueva sin mutantes.
        Si la matriz ya está sana, no realiza cambios.

        Args:
            matriz (List[List[str]]): Matriz de ADN a sanar.

        Returns:
            Tuple[List[List[str]], bool]: Matriz de ADN (nueva o la misma) y un booleano indicando si se realizaron cambios.
        """
        detector = Detector(matriz)
        if not detector.detectar_mutantes():
            return matriz, False  # La matriz ya está sana

        try:
            nueva_matriz = [
                [random.choice("ATCG") for _ in range(len(matriz[0]))]
                for _ in range(len(matriz))
            ]
            while Detector(nueva_matriz).detectar_mutantes():
                nueva_matriz = [
                    [random.choice("ATCG") for _ in range(len(matriz[0]))]
                    for _ in range(len(matriz))
                ]
            return nueva_matriz, True  # La matriz fue sanada
        except Exception as e:
            print(f"Error al sanar la matriz: {e}")
            return matriz, False
