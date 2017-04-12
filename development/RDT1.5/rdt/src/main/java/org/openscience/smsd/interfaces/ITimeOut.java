/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.smsd.interfaces;

/**
 *
 * @author Asad
 */
public interface ITimeOut {

    /**
     * get timeout in mins for bond insensitive searches
     * @return the bondInSensitive TimeOut
     */
    double getBondInSensitiveCDKMCSTimeOut();

    /**
     * get timeout in mins for bond insensitive searches
     * @return the bondInSensitive TimeOut
     */
    double getBondInSensitiveMCSPlusTimeOut();

    /**
     * get timeout in mins for bond insensitive searches
     * @return the bondInSensitive TimeOut
     */
    double getBondInSensitiveVFTimeOut();

    /**
     * get timeout in mins for bond sensitive searches
     * @return the bondSensitive TimeOut
     */
    double getBondSensitiveCDKMCSTimeOut();

    /**
     * get timeout in mins for bond sensitive searches
     * @return the bondSensitive TimeOut
     */
    double getBondSensitiveMCSPlusTimeOut();

    /**
     * get timeout in mins for bond sensitive searches
     * @return the bondSensitive TimeOut
     */
    double getBondSensitiveVFTimeOut();

    /**
     * set timeout in mins (default 1.00 min) for bond insensitive searches
     * @param bondInSensitiveTimeOut the bond insensitive
     */
    void setBondInSensitiveCDKMCSTimeOut(double bondInSensitiveTimeOut);

    /**
     * set timeout in mins (default 1.00 min) for bond insensitive searches
     * @param bondInSensitiveTimeOut the bond insensitive
     */
    void setBondInSensitiveMCSPlusTimeOut(double bondInSensitiveTimeOut);

    /**
     * set timeout in mins (default 1.00 min) for bond insensitive searches
     * @param bondInSensitiveTimeOut the bond insensitive
     */
    void setBondInSensitiveVFTimeOut(double bondInSensitiveTimeOut);

    /**
     * set timeout in mins (default 0.10 min) for bond sensitive searches
     * @param bondSensitiveTimeOut the bond Sensitive Timeout in mins (default 0.30 min)
     */
    void setBondSensitiveCDKMCSTimeOut(double bondSensitiveTimeOut);

    /**
     * set timeout in mins (default 0.10 min) for bond sensitive searches
     * @param bondSensitiveTimeOut the bond Sensitive Timeout in mins (default 0.30 min)
     */
    void setBondSensitiveMCSPlusTimeOut(double bondSensitiveTimeOut);

    /**
     * set timeout in mins (default 0.10 min) for bond sensitive searches
     * @param bondSensitiveTimeOut the bond Sensitive Timeout in mins (default 0.30 min)
     */
    void setBondSensitiveVFTimeOut(double bondSensitiveTimeOut);

}
