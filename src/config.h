/** @file config.h
 * @brief Header file for compile time settings
 *
 * @warning Recompile entire code whenever you make any changes to this file:
 @code
   make clean; make
 @endcode
 */

/* Use color highlighting of (some) output */
#define USECOLOR 1

/* Edit the following to save (some_ paths generated by the
   single level method to disk */

/* Save (some) states to disk? */
//#define SAVE_STATES 1

/* First state to save */
#define SAVE_FIRST_STATE 5000

/* Last path to save */
#define SAVE_LAST_STATE 5100

/* Log QoI in singlelevel method? */
//#define LOG_QOI 1
