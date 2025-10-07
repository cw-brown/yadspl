/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a struct for SCHP5 packets
 * @version 0.1
 * @date 2025-10-02
 */
# ifndef PACKET_H
# define PACKET_H

#include <cstdint>
#include <vector>

/**
 * @brief Implements a packet compliant to the SCHP5 protocol
 */
struct Packet {

    // Fields
    /**
     * @brief The flag at the start and stop of each packet
     */
    static const uint8_t flag = 0b01111110;

    /**
     * @brief Describes the packet type
     */
    uint8_t controlCode;

    /**
     * @brief Describes how many bytes are contained in the packet
     */
    uint8_t dataLength;

    /**
     * @brief The address of the sender of this packet
     */
    uint32_t recieverAddress;

    /**
     * @brief The address of the sender of this packet
     */
    uint32_t senderAddress;

    /**
     * @brief The sequence number of the current packet
     */
    uint8_t sequenceNum;

    /** 
     * @brief The data content of this packet
    */
   bool * data;

   /**
    * @brief Contains ISO-3309 CRC-32 checksum
    */
   uint32_t erc;

};

#endif