# Библиотека OpenVX Ext

Тестовое задание для студентов: реализация функций видеоаналитики + демонстрационные примеры.

Документ с заданием находится в [Description/README.md](Description/README.md).

# Структура библиотеки

## Обозначения

__[Название папки]__ - обозначает каталог

_Название файла_ - обозначает файл

## openvx_ext

* __[Image]__ - Содержит вспомогательные изображения
* __[Lib]__ - Содержит исходный код алгоритмов и вспомогательных функций
    * __[Common]__ - Содержит вспомогательные функции для реализации алгоритмов, заголовочные файлы с типами данных и функции runtime
        * __[openvx]__ - Содержит заголовочные файлы OpenVX
            * _vx.h_ - Заголовочный файл OpenVX
            * _vx\_api.h_ - Определение API для OpenVX
            * ...
        * ...
    * __[Kernels]__ - Содержит исходный код различных реализаций алгоритма
        * __[ref]__ - Содержит исходный код эталонных функций под RISC
            * _BoxFilterBordersTemplate.h_ - Содержит функции обработки границ для BoxFilte
            * _ref\_AbsDiff.c_ - Содержит эталонную реализацию функции Absolute Difference
            * ...
        * _ref.h_ - Заголовочный файл с объявлениями эталонных функций под RISC
* __[Demo]__ - Демонстрационный проект
