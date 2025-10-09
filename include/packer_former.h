/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for forming SCHP5 packets
 * @version 0.1
 * @date 2025-10-02
 */
# ifndef PACKET_FORMER_H
# define PACKET_FORMER_H

#include "packet.h"
#include "ring.hpp"
#include "crcinterface.h"
#include <cstdint>
#include <queue>

/**
 * @brief Forms packets according to the SCHP5 protocol for packet radios
 */
class PacketFormer {

    private:

    /**
    * @brief The address of the sender
    */
    const uint32_t senderAddress;

    /**
    * @brief The address of the reciever
    */
    const uint32_t recieverAddress;

    /**
    * @brief The current length of the data portion of the packets
    */
    uint8_t dataLength;

    /**
     * @brief The number of the current packet in sequence
     */
    uint8_t sequenceNumber;

    crcutil_interface::CRC * crcGenny;

    public:

    /**
     * @brief Pointer to the buffer where the packet former reads in arbitrary binary data formatted as a buffer of bytes
     */
    //const std::ring<uint8_t> * inputBuffer;

    /**
     * @brief Pointer to the buffer where the formed packets are stored
     */
    //const std::ring<Packet> * outputBuffer;

    // For testing purposes due to issues with ring class

     /**
     * @brief Pointer to the buffer where the packet former reads in arbitrary binary data formatted as a buffer of bytes
     */
    std::queue<uint8_t> * inputBuffer;

    /**
     * @brief Pointer to the buffer where the formed packets are stored
     */
    std::queue<Packet> * outputBuffer;

    public:

    /**
     * @brief Construct a new packet former
     */
    PacketFormer(
        std::ring<uint8_t> * input_buffer,
        std::ring<Packet> * output_buffer,
        uint32_t sender_address,
        uint32_t reciever_address,
        uint8_t data_length);


    /**
     * @brief Set Data Length
     */
    void setDataLength(uint8_t data_length);

    /**
     * @brief Set Data Length
     */
    void resetSequenceNumber();

    /**
     * @brief Form the next packet
     */
    void formNextPacket();

};

#endif